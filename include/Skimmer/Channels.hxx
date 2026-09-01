#pragma once

#include <cstddef>
#include <span>

#include <Math/Point3D.h>

#include "common/Constants.hpp"
#include "common/Schema_FoundHdibaryon.hpp"
#include "common/Schema_FoundSexaquark.hpp"

#include "Skimmer/Process.hxx"
#include "Skimmer/VariableRegistry.hxx"

namespace Skimmer {

// ## Channel Traits ## //

// One struct per analysis channel, giving the generic skim everything it needs to know: which input
// schema to read, how to turn index `i` of the index-parallel candidate vectors into a `Cached::`
// object, how to label that candidate, and how many signals were injected per event.
//
// The RNTuple name is *not* here: it comes from the config's `samples[].ntuple`, so that the config
// stays the one place a file and what to read out of it are named together.
//
// `Size` and `McSize` exist as a pair because `Label` indexes the MC vector with an index bounded by
// the reconstructed one. That parallelism is a `t2ds` invariant rather than anything this program can
// enforce, so the skim compares the two once per event: a mismatch is an out-of-bounds read, and a
// silent one is far worse than a stopped job.

struct TraitsChannelA {
    using Schema = ::Schema::FoundSexaquark;
    using Cached = ::Cached::ChannelA;

    static constexpr std::span<const Variables::Definition<Cached>> kVariables{Variables::DB_ChannelA};
    static constexpr unsigned int kNInjectedPerEvent{E2T::NSexaReactionsPerEvent};

    static std::size_t Size(const Schema& schema) { return schema.ChannelA.size(); }
    static std::size_t McSize(const Schema& schema) { return schema.MC_ChannelA.size(); }

    static Cached Build(const Schema& schema, std::size_t i, const ROOT::Math::XYZPoint& pv) {
        return Cached{schema.ChannelA[i], schema.ChannelA_V0A[i], schema.ChannelA_V0B[i], pv};
    }

    static Process::EProcess Label(const Schema& schema, std::size_t i, const Cached& cached) {
        if (cached.IsWrongSignChannel) return Process::kReference;
        const auto& mc = schema.MC_ChannelA[i];
        return Process::Classify(mc.IsTrueSignal, mc.IsRealBkg);
    }
};

struct TraitsChannelD {
    using Schema = ::Schema::FoundSexaquark;
    using Cached = ::Cached::ChannelD;

    static constexpr std::span<const Variables::Definition<Cached>> kVariables{Variables::DB_ChannelD};
    static constexpr unsigned int kNInjectedPerEvent{E2T::NSexaReactionsPerEvent};

    static std::size_t Size(const Schema& schema) { return schema.ChannelD.size(); }
    static std::size_t McSize(const Schema& schema) { return schema.MC_ChannelD.size(); }

    static Cached Build(const Schema& schema, std::size_t i, const ROOT::Math::XYZPoint& pv) {
        return Cached{schema.ChannelD[i], schema.ChannelD_V0[i], pv};
    }

    static Process::EProcess Label(const Schema& schema, std::size_t i, const Cached& cached) {
        if (cached.IsWrongSignChannel) return Process::kReference;
        const auto& mc = schema.MC_ChannelD[i];
        return Process::Classify(mc.IsTrueSignal, mc.IsRealBkg);
    }
};

struct TraitsChannelH {
    using Schema = ::Schema::FoundSexaquark;
    using Cached = ::Cached::ChannelH;

    static constexpr std::span<const Variables::Definition<Cached>> kVariables{Variables::DB_ChannelH};
    static constexpr unsigned int kNInjectedPerEvent{E2T::NSexaReactionsPerEvent};

    static std::size_t Size(const Schema& schema) { return schema.ChannelH.size(); }
    static std::size_t McSize(const Schema& schema) { return schema.MC_ChannelH.size(); }

    static Cached Build(const Schema& schema, std::size_t i, const ROOT::Math::XYZPoint& pv) { return Cached{schema.ChannelH[i], pv}; }

    static Process::EProcess Label(const Schema& schema, std::size_t i, const Cached& cached) {
        if (cached.IsWrongSignChannel) return Process::kReference;
        const auto& mc = schema.MC_ChannelH[i];
        return Process::Classify(mc.IsTrueSignal, mc.IsRealBkg);
    }
};

struct TraitsHdibaryon {
    using Schema = ::Schema::FoundHdibaryon;
    using Cached = ::Cached::Hdibaryon;

    static constexpr std::span<const Variables::Definition<Cached>> kVariables{Variables::DB_Hdibaryon};
    static constexpr unsigned int kNInjectedPerEvent{E2T::NInjectedHdibaryonsPerEvent};

    static std::size_t Size(const Schema& schema) { return schema.Hdibaryon.size(); }
    static std::size_t McSize(const Schema& schema) { return schema.MC_Hdibaryon.size(); }

    static Cached Build(const Schema& schema, std::size_t i, const ROOT::Math::XYZPoint& pv) {
        return Cached{schema.Hdibaryon[i], schema.Lambda1[i], schema.Lambda2[i], pv};
    }

    static Process::EProcess Label(const Schema& schema, std::size_t i, const Cached& cached) {
        if (cached.IsMixedChannel()) return Process::kReference;
        const auto& mc = schema.MC_Hdibaryon[i];
        return Process::Classify(mc.IsTrueSignal, mc.IsRealBkg);
    }
};

}  // namespace Skimmer
