#ifndef T2S_CONSTANTS_HXX
#define T2S_CONSTANTS_HXX

#include <string_view>

namespace Tree2Secondaries {

enum class ReactionChannel : char { A = 'A', D = 'D', E = 'E', H = 'H', AntiA = 'a', AntiD = 'd', AntiE = 'e', AntiH = 'h', All = ' ' };

namespace Name {
static constexpr std::string_view Events{"Events"};
static constexpr std::string_view PackedEvents{"PackedEvents"};
static constexpr std::string_view ChannelA{"ChannelA"};
static constexpr std::string_view AntiA{"AntiA"};
static constexpr std::string_view ChannelD{"ChannelD"};
static constexpr std::string_view AntiD{"AntiD"};
static constexpr std::string_view ChannelE{"ChannelE"};
static constexpr std::string_view AntiE{"AntiE"};
static constexpr std::string_view ChannelH{"ChannelH"};
static constexpr std::string_view AntiH{"AntiH"};
}  // namespace Name

namespace Acronym {
static constexpr std::string_view AntiLambda{"AL"};
static constexpr std::string_view Lambda{"L"};
static constexpr std::string_view KaonZeroShort{"K0S"};
static constexpr std::string_view AntiProton{"AP"};
static constexpr std::string_view Proton{"P"};
static constexpr std::string_view NegKaon{"NK"};
static constexpr std::string_view PosKaon{"PK"};
static constexpr std::string_view PiMinus{"PM"};
static constexpr std::string_view PiPlus{"PP"};
}  // namespace Acronym

enum class PdgCode : short {
    AntiLambda = -3122,
    Lambda = 3122,
    KaonZeroShort = 310,
    AntiNeutron = -2112,
    Neutron = 2112,
    AntiProton = -2212,
    Proton = 2212,
    NegKaon = -321,
    PosKaon = 321,
    PiMinus = -211,
    PiPlus = 211
};

namespace PdgMass {
static constexpr double Lambda{1.1156830};          // (GeV/c^2)
static constexpr double KaonZeroShort{0.49761100};  // (GeV/c^2)
static constexpr double Neutron{0.93956540};        // (GeV/c^2)
static constexpr double Proton{0.93827210};         // (GeV/c^2)
static constexpr double Kaon{0.49367700};           // (GeV/c^2)
static constexpr double Pion{0.13957040};           // (GeV/c^2)
}  // namespace PdgMass

namespace Const {
static constexpr double Kappa{0.000299792458};  // (GeV/c) / (kG/cm)
static constexpr double AbsAlmostZero{1.E-8};
static constexpr double BigNumber{1.E8};
static constexpr int DummyInt{-1};
static constexpr float DummyFloat{-999.};
static constexpr double DummyDouble{-999.};
}  // namespace Const

}  // namespace Tree2Secondaries

#endif  // T2S_CONSTANTS_HXX
