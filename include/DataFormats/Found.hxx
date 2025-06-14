#ifndef T2S_STRUCTS_FOUND_HXX
#define T2S_STRUCTS_FOUND_HXX

#include <cstddef>

namespace Tree2Secondaries::Found {

struct alignas(32) State {
    size_t Entry{};
    float X{};
    float Y{};
    float Z{};
    float Px{};
    float Py{};
    float Pz{};
    float E{};
};

struct alignas(32) Sexaquark : State {
    unsigned int RunNumber{};
    unsigned int DirNumber{};
    unsigned int DirNumberB{};
    unsigned int EventNumber{};
    float E_MinusNucleon{};
};

struct alignas(32) ChannelA : Sexaquark {
    State V0A;
    State V0B;

    size_t V0A_Neg_Entry{};
    size_t V0A_Pos_Entry{};
    size_t V0B_Neg_Entry{};
    size_t V0B_Pos_Entry{};

    float V0A_DecayX{};
    float V0A_DecayY{};
    float V0A_DecayZ{};

    float V0B_DecayX{};
    float V0B_DecayY{};
    float V0B_DecayZ{};
};

struct alignas(32) ChannelD : Sexaquark {
    State V0;
    State Kaon;

    size_t V0_Neg_Entry{};
    size_t V0_Pos_Entry{};

    float V0_DecayX{};
    float V0_DecayY{};
    float V0_DecayZ{};
};

struct alignas(32) ChannelE : ChannelD {
    State PiMinus;
    State PiPlus;
};

struct alignas(32) ChannelH : Sexaquark {
    State Kaon1;
    State Kaon2;
};

}  // namespace Tree2Secondaries::Found

#endif  // T2S_STRUCTS_FOUND_HXX
