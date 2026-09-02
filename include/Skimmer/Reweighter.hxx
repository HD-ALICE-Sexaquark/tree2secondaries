#pragma once

#include <cstddef>
#include <cstdint>
#include <random>
#include <string_view>
#include <vector>

#include <TH1.h>

namespace Skimmer {

// ## The Shape Weight ## //

// The dedicated antisexaquark MC injects flat in pt over [0, 5) GeV/c and flat in 2D radius over [5, 180) cm (`common/docs/MC_PRODUCTIONS.md`,
// `E2T::InjectedAntiSexa_*`, `common/Constants.hpp`). Neither is physical, so a signal candidate has to be reweighted into the shape it would have
// had:
//
//   pt      blast-wave templates, one per centrality class, keyed "<mass>_<class>" (e.g. "1.80_10-20")
//   radius  one "Radius2D" template, the real-data material budget measured with photon conversions
//
// Each factor is `template(x) / flat-prior-density`, restricted to the injected range, so the mean weight over that range is 1. That makes this a
// shape weight and nothing else: it says what the distribution should look like, not how many antisexaquarks a real run would hold.
//
// The yield weight should be applied by the consumer of the cache.

class Reweighter {
   public:
    Reweighter(const Reweighter&) = delete;
    Reweighter(Reweighter&&) = delete;
    Reweighter& operator=(const Reweighter&) = delete;
    Reweighter& operator=(Reweighter&&) = delete;
    ~Reweighter() = default;

    // Throws if a named file cannot be opened or does not hold every template `signal_mass` asks for.
    Reweighter(double signal_mass, std::string_view pt_path, std::string_view radius_path);

    // `pt` and `radius` are MC truth, taken from the injected antisexaquark linked to the candidate.
    // NOTE: not const -- an event whose centrality falls outside the classes draws a random one, and both
    //       that draw and its tally are state.
    [[nodiscard]] double Weight(float centrality, double pt, double radius);

    [[nodiscard]] bool IsActive() const { return fUsePt || fUseRadius; }
    [[nodiscard]] std::uint64_t NRandomizedCentrality() const { return fNRandomizedCentrality; }

   private:
    void LoadPt(double signal_mass, std::string_view path);
    void LoadRadius(std::string_view path);

    // Maps `centrality` onto a class, drawing a random one when it falls outside [0, 90).
    [[nodiscard]] std::size_t ResolveCentralityClass(float centrality);

    [[nodiscard]] double WeightPt(std::size_t cc, double pt) const;
    [[nodiscard]] double WeightRadius(double radius) const;

    // -- pt
    std::vector<TH1D> fTemplatePt;  // one per centrality class, index-parallel to `DB::Centrality::Name`
    std::vector<double> fScalePt;   // parallel to `fTemplatePt`
    bool fUsePt{false};
    // -- radius
    TH1D fTemplateRadius;  // single histogram
    double fScaleRadius{1.};
    bool fUseRadius{false};
    // -- the centrality fallback
    std::mt19937 fRandom{0};  // a fixed seed, so that re-running a skim reproduces its cache.
    std::uint64_t fNRandomizedCentrality{0};
};

}  // namespace Skimmer
