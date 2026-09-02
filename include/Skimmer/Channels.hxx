#pragma once

#include <cstddef>
#include <span>
#include <string_view>

#include <Math/Point3D.h>

#include "common/Cached_InjectedSexa.hpp"
#include "common/Constants.hpp"
#include "common/Schema_FoundHdibaryon.hpp"
#include "common/Schema_FoundSexaquark.hpp"

#include "Skimmer/Classification.hxx"
#include "Skimmer/VariableRegistry.hxx"

namespace Skimmer {

// ## Channel Traits ## //

// One struct per analysis channel, giving the generic skim everything it needs to know: which input
// schema to read, which RNTuple that schema lives in and which one it caches into, how to turn index `i`
// of the index-parallel candidate vectors into a `Cached::` object, and how to label that candidate.
//
// `Size` and `McSize` exist as a pair because `Label` indexes the MC vector with an index bounded by the
// reconstructed one. That parallelism is a `t2ds` invariant rather than anything this program can enforce,
// so the skim compares the two once per event: a mismatch is an out-of-bounds read, and a silent one is
// far worse than a stopped job.

struct TraitsChannelA {
    using Schema = ::Schema::FoundSexaquark;
    using Cached = ::Cached::ChannelA;

    static constexpr std::span<const Variables::Definition<Cached>> kVariables{Variables::DB_ChannelA};
    static constexpr std::string_view kName_InputRNT{T2DS::Name_FoundSexaquarkRNT};
    static constexpr std::string_view kName_OutputRNT{Skimmer::Name_CachedSexaquarkRNT};

    static constexpr bool kHasInjectedSexa = true;  // needed for reweighting

    static std::size_t Size(const Schema& schema) { return schema.ChannelA.size(); }
    static std::size_t McSize(const Schema& schema) { return schema.MC_ChannelA.size(); }

    static Cached Build(const Schema& schema, std::size_t i, const ROOT::Math::XYZPoint& pv) {
        return Cached{schema.ChannelA[i], schema.ChannelA_V0A[i], schema.ChannelA_V0B[i], pv};
    }

    // The MC PV, not the reconstructed one: this is an entirely generator-level object, and only the truth
    // vertex belongs in it -- even though neither `Pt()` nor `SV_Radius2D()` reads the reference point.
    static ::Cached::InjectedSexa Injected(const Schema& schema, std::size_t i) {
        return {schema.MC_ChannelA[i], {schema.MC_Event.PV_X, schema.MC_Event.PV_Y, schema.MC_Event.PV_Z}};
    }

    static Classification::EClassification Label(const Schema& schema, std::size_t i, const Cached& cached) {
        const auto& mc = schema.MC_ChannelA[i];
        if (cached.IsWrongSignChannel) return mc.IsRealBkg ? Classification::kReferenceReal : Classification::kReferenceHybrid;
        return Classification::Classify(mc.IsTrueSignal, mc.IsRealBkg);
    }
};

struct TraitsChannelD {
    using Schema = ::Schema::FoundSexaquark;
    using Cached = ::Cached::ChannelD;

    static constexpr std::span<const Variables::Definition<Cached>> kVariables{Variables::DB_ChannelD};
    static constexpr std::string_view kName_InputRNT{T2DS::Name_FoundSexaquarkRNT};
    static constexpr std::string_view kName_OutputRNT{Skimmer::Name_CachedSexaquarkRNT};

    static constexpr bool kHasInjectedSexa = true;  // needed for reweighting

    static std::size_t Size(const Schema& schema) { return schema.ChannelD.size(); }
    static std::size_t McSize(const Schema& schema) { return schema.MC_ChannelD.size(); }

    static Cached Build(const Schema& schema, std::size_t i, const ROOT::Math::XYZPoint& pv) {
        return Cached{schema.ChannelD[i], schema.ChannelD_V0[i], pv};
    }

    static ::Cached::InjectedSexa Injected(const Schema& schema, std::size_t i) {
        return {schema.MC_ChannelD[i], {schema.MC_Event.PV_X, schema.MC_Event.PV_Y, schema.MC_Event.PV_Z}};
    }

    static Classification::EClassification Label(const Schema& schema, std::size_t i, const Cached& cached) {
        const auto& mc = schema.MC_ChannelD[i];
        if (cached.IsWrongSignChannel) return mc.IsRealBkg ? Classification::kReferenceReal : Classification::kReferenceHybrid;
        return Classification::Classify(mc.IsTrueSignal, mc.IsRealBkg);
    }
};

struct TraitsChannelH {
    using Schema = ::Schema::FoundSexaquark;
    using Cached = ::Cached::ChannelH;

    static constexpr std::span<const Variables::Definition<Cached>> kVariables{Variables::DB_ChannelH};
    static constexpr std::string_view kName_InputRNT{T2DS::Name_FoundSexaquarkRNT};
    static constexpr std::string_view kName_OutputRNT{Skimmer::Name_CachedSexaquarkRNT};

    static constexpr bool kHasInjectedSexa = true;  // needed for reweighting

    static std::size_t Size(const Schema& schema) { return schema.ChannelH.size(); }
    static std::size_t McSize(const Schema& schema) { return schema.MC_ChannelH.size(); }

    static Cached Build(const Schema& schema, std::size_t i, const ROOT::Math::XYZPoint& pv) { return Cached{schema.ChannelH[i], pv}; }

    static ::Cached::InjectedSexa Injected(const Schema& schema, std::size_t i) {
        return {schema.MC_ChannelH[i], {schema.MC_Event.PV_X, schema.MC_Event.PV_Y, schema.MC_Event.PV_Z}};
    }

    static Classification::EClassification Label(const Schema& schema, std::size_t i, const Cached& cached) {
        const auto& mc = schema.MC_ChannelH[i];
        if (cached.IsWrongSignChannel) return mc.IsRealBkg ? Classification::kReferenceReal : Classification::kReferenceHybrid;
        return Classification::Classify(mc.IsTrueSignal, mc.IsRealBkg);
    }
};

struct TraitsLambdaPair {
    using Schema = ::Schema::FoundHdibaryon;
    using Cached = ::Cached::Hdibaryon;

    static constexpr std::span<const Variables::Definition<Cached>> kVariables{Variables::DB_Hdibaryon};
    static constexpr std::string_view kName_InputRNT{T2DS::Name_FoundHdibaryonRNT};
    static constexpr std::string_view kName_OutputRNT{Skimmer::Name_CachedHdibaryonRNT};

    static constexpr bool kHasInjectedSexa = false;  // needed for reweighting

    static std::size_t Size(const Schema& schema) { return schema.Hdibaryon.size(); }
    static std::size_t McSize(const Schema& schema) { return schema.MC_Hdibaryon.size(); }

    static Cached Build(const Schema& schema, std::size_t i, const ROOT::Math::XYZPoint& pv) {
        return Cached{schema.Hdibaryon[i], schema.Lambda1[i], schema.Lambda2[i], pv};
    }

    static Classification::EClassification Label(const Schema& schema, std::size_t i, const Cached& cached) {
        const auto& mc = schema.MC_Hdibaryon[i];
        if (cached.IsMixedChannel()) return mc.IsRealBkg ? Classification::kReferenceReal : Classification::kReferenceHybrid;
        return Classification::Classify(mc.IsTrueSignal, mc.IsRealBkg);
    }
};

}  // namespace Skimmer
