#include <cstddef>
#include <format>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <string_view>

#include <TFile.h>
#include <TH1.h>

#include "common/Constants.hpp"
#include "common/DB_Centrality.hpp"

#include "Skimmer/Reweighter.hxx"

namespace Skimmer {

// # Construction # //

Reweighter::Reweighter(double signal_mass, std::string_view pt_path, std::string_view radius_path)
    : fUsePt{!pt_path.empty()}, fUseRadius{!radius_path.empty()} {
    LoadPt(signal_mass, pt_path);
    LoadRadius(radius_path);
}

void Reweighter::LoadPt(double signal_mass, std::string_view path) {

    if (!fUsePt) return;

    const std::unique_ptr<TFile> file{TFile::Open(std::string(path).c_str(), "READ")};
    if (!file || file->IsZombie()) {
        throw std::runtime_error(std::format("could not open pt weights file \"{}\"", path));
    }

    fTemplatePt.reserve(DB::Centrality::NClasses);
    for (std::size_t cc = 0; cc < DB::Centrality::NClasses; ++cc) {
        const std::string name = std::format("{:.2f}_{}", signal_mass, DB::Centrality::Name[cc]);
        const auto* histogram = file->Get<TH1D>(name.c_str());
        if (!histogram) {
            throw std::runtime_error(std::format("\"{}\" has no TH1D \"{}\" -- it holds no blast-wave template for signal mass {:.3f}",  //
                                                 path, name, signal_mass));
        }
        fTemplatePt.emplace_back(*histogram);
        fTemplatePt.back().SetDirectory(nullptr);
    }

    fScalePt.reserve(fTemplatePt.size());
    for (const auto& tmpl : fTemplatePt) {
        const double integral = tmpl.Integral(1, tmpl.FindFixBin(E2T::InjectedAntiSexa_MaxPt) - 1, "width");
        if (integral <= 0.) {
            throw std::runtime_error(std::format("\"{}\": template \"{}\" is empty over the injected pt range [{}, {}) GeV/c",  //
                                                 path, tmpl.GetName(), E2T::InjectedAntiSexa_MinPt, E2T::InjectedAntiSexa_MaxPt));
        }
        fScalePt.push_back((E2T::InjectedAntiSexa_MaxPt - E2T::InjectedAntiSexa_MinPt) / integral);
    }
}

void Reweighter::LoadRadius(std::string_view path) {

    if (!fUseRadius) return;

    const std::unique_ptr<TFile> file{TFile::Open(std::string(path).c_str(), "READ")};
    if (!file || file->IsZombie()) {
        throw std::runtime_error(std::format("could not open radius weights file \"{}\"", path));
    }

    const auto* histogram = file->Get<TH1D>("Radius2D");
    if (!histogram) {
        throw std::runtime_error(std::format("\"{}\" has no TH1D \"Radius2D\" -- it holds no material budget template", path));
    }
    fTemplateRadius = *histogram;
    fTemplateRadius.SetDirectory(nullptr);

    const double integral = fTemplateRadius.Integral(fTemplateRadius.FindFixBin(E2T::InjectedAntiSexa_MinRadius),  //
                                                     fTemplateRadius.GetNbinsX(), "width");
    if (integral <= 0.) {
        throw std::runtime_error(std::format("\"{}\": template \"Radius2D\" is empty over the injected radius range [{}, {}) cm",  //
                                             path, E2T::InjectedAntiSexa_MinRadius, E2T::InjectedAntiSexa_MaxRadius));
    }
    fScaleRadius = (E2T::InjectedAntiSexa_MaxRadius - E2T::InjectedAntiSexa_MinRadius) / integral;
}

// # Lookup # //

std::size_t Reweighter::ResolveCentralityClass(float centrality) {

    const std::size_t cc = DB::Centrality::ClassOf(centrality);
    if (cc < DB::Centrality::NClasses) return cc;

    // -- if centrality out-of-range, get random centrality
    ++fNRandomizedCentrality;
    std::uniform_real_distribution<float> flat{DB::Centrality::Edges.front(), DB::Centrality::Edges.back()};
    const std::size_t drawn = DB::Centrality::ClassOf(flat(fRandom));
    return drawn < DB::Centrality::NClasses ? drawn : DB::Centrality::NClasses - 1;
}

double Reweighter::WeightPt(std::size_t cc, double pt) const {

    // loudly reject out-of-range truth pt
    if (pt < E2T::InjectedAntiSexa_MinPt || pt >= E2T::InjectedAntiSexa_MaxPt) {
        throw std::runtime_error(std::format("injected pt {:.4f} GeV/c lies outside the injected range [{}, {}) GeV/c",  //
                                             pt, E2T::InjectedAntiSexa_MinPt, E2T::InjectedAntiSexa_MaxPt));
    }

    const TH1D& tmpl = fTemplatePt[cc];
    return fScalePt[cc] * tmpl.GetBinContent(tmpl.FindFixBin(pt));
}

double Reweighter::WeightRadius(double radius) const {

    // loudly reject out-of-range truth 2d radius
    if (radius < E2T::InjectedAntiSexa_MinRadius || radius >= E2T::InjectedAntiSexa_MaxRadius) {
        throw std::runtime_error(std::format("injected 2D radius {:.4f} cm lies outside the injected range [{}, {}) cm",  //
                                             radius, E2T::InjectedAntiSexa_MinRadius, E2T::InjectedAntiSexa_MaxRadius));
    }

    return fScaleRadius * fTemplateRadius.GetBinContent(fTemplateRadius.FindFixBin(radius));
}

double Reweighter::Weight(float centrality, double pt, double radius) {
    const double weight_pt = fUsePt ? WeightPt(ResolveCentralityClass(centrality), pt) : 1.;
    const double weight_radius = fUseRadius ? WeightRadius(radius) : 1.;
    return weight_pt * weight_radius;
}

}  // namespace Skimmer
