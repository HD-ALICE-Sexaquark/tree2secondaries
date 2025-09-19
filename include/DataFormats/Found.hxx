#ifndef T2S_STRUCTS_FOUND_HXX
#define T2S_STRUCTS_FOUND_HXX

#include "Math/Constants.hxx"

namespace Tree2Secondaries::Found {

struct alignas(T2S_SIMD_ALIGN) State {
    float X{};
    float Y{};
    float Z{};
    float Px{};
    float Py{};
    float Pz{};
    float E{};
};

struct alignas(T2S_SIMD_ALIGN) Sexaquark : State {
    unsigned int RunNumber{};
    unsigned int DirNumber{};
    unsigned int DirNumberB{};
    unsigned int EventNumber{};
    float MagneticField{};
    float PV_Xv{};
    float PV_Yv{};
    float PV_Zv{};
    float MC_PV_Xv{};
    float MC_PV_Yv{};
    float MC_PV_Zv{};
    float E_MinusNucleon{};
};

struct alignas(T2S_SIMD_ALIGN) ChannelA : Sexaquark {
    State V0A;
    State V0B;

    int V0A_Entry{};
    int V0A_Neg_Entry{};
    int V0A_Pos_Entry{};

    float V0A_X_AtPCA{};
    float V0A_Y_AtPCA{};
    float V0A_Z_AtPCA{};

    int V0B_Entry{};
    int V0B_Neg_Entry{};
    int V0B_Pos_Entry{};

    float V0B_X_AtPCA{};
    float V0B_Y_AtPCA{};
    float V0B_Z_AtPCA{};
};

struct alignas(T2S_SIMD_ALIGN) ChannelD : Sexaquark {
    State V0;
    State Kaon;

    int V0_Entry{};
    int V0_Neg_Entry{};
    int V0_Pos_Entry{};

    float V0_X_AtPCA{};
    float V0_Y_AtPCA{};
    float V0_Z_AtPCA{};

    int Kaon_Entry{};

    float Kaon_X_AtPCA{};
    float Kaon_Y_AtPCA{};
    float Kaon_Z_AtPCA{};
    float Kaon_Px_AtPCA{};
    float Kaon_Py_AtPCA{};
    float Kaon_Pz_AtPCA{};
};

struct alignas(T2S_SIMD_ALIGN) ChannelE : ChannelD {
    State PiMinus;
    State PiPlus;
};

struct alignas(T2S_SIMD_ALIGN) ChannelH : Sexaquark {
    State Kaon1;
    State Kaon2;
};

struct alignas(T2S_SIMD_ALIGN) Injected : State {
    unsigned int RunNumber{};
    unsigned int DirNumber{};
    unsigned int EventNumber{};
    int ReactionID{};
};

struct alignas(T2S_SIMD_ALIGN) MC_Sexaquark {
    float X{};
    float Y{};
    float Z{};
    float BeforePx{};
    float BeforePy{};
    float BeforePz{};
    float BeforeE{};

    float NucleonPx{};
    float NucleonPy{};
    float NucleonPz{};
    float NucleonE{};

    float AfterPx{};
    float AfterPy{};
    float AfterPz{};
    float AfterE{};

    int ReactionID{};
    bool IsSignal{};
    bool IsHybrid{};
};

struct alignas(T2S_SIMD_ALIGN) MC_ChannelA : MC_Sexaquark {
    float V0A_Px{};
    float V0A_Py{};
    float V0A_Pz{};
    float V0A_E{};

    float V0B_Px{};
    float V0B_Py{};
    float V0B_Pz{};
    float V0B_E{};

    int V0A_Entry{};
    int V0A_Mother_Entry{};
    int V0A_Neg_Entry{};
    int V0A_Pos_Entry{};
    int V0A_PdgCode{};
    int V0A_Mother_PdgCode{};
    int V0A_ReactionID{};

    int V0B_Entry{};
    int V0B_Mother_Entry{};
    int V0B_Neg_Entry{};
    int V0B_Pos_Entry{};
    int V0B_PdgCode{};
    int V0B_Mother_PdgCode{};
    int V0B_ReactionID{};

    bool V0A_IsTrue{};
    bool V0A_IsSignal{};
    bool V0A_IsSecondary{};
    bool V0A_IsHybrid{};

    bool V0B_IsTrue{};
    bool V0B_IsSignal{};
    bool V0B_IsSecondary{};
    bool V0B_IsHybrid{};
};

struct alignas(T2S_SIMD_ALIGN) MC_ChannelD : MC_Sexaquark {
    float V0_Px{};
    float V0_Py{};
    float V0_Pz{};
    float V0_E{};

    float Kaon_Px{};
    float Kaon_Py{};
    float Kaon_Pz{};
    float Kaon_E{};

    int V0_Entry{};
    int V0_Neg_Entry{};
    int V0_Pos_Entry{};
    int V0_PdgCode{};
    int V0_Mother_Entry{};
    int V0_Mother_PdgCode{};
    int V0_ReactionID{};

    int Kaon_Entry{};
    int Kaon_PdgCode{};
    int Kaon_Mother_Entry{};
    int Kaon_Mother_PdgCode{};
    int Kaon_GrandMother_Entry{};
    int Kaon_GrandMother_PdgCode{};
    int Kaon_ReactionID{};

    bool V0_IsTrue{};
    bool V0_IsSignal{};
    bool V0_IsSecondary{};
    bool V0_IsHybrid{};

    bool Kaon_IsTrue{};
    bool Kaon_IsSignal{};
    bool Kaon_IsSecondary{};
    bool Kaon_IsHybrid{};
};

}  // namespace Tree2Secondaries::Found

#endif  // T2S_STRUCTS_FOUND_HXX
