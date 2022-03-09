#ifndef OPTIONS_H
#define OPTIONS_H

namespace multipole_conv {
enum class MPOptions {
  none = 1 << 0,
  complex = 1 << 1,
  normalisation = 1 << 2,
  condon_shortley_phase = 1 << 3
};

constexpr MPOptions operator|(MPOptions f1, MPOptions f2) {
  return static_cast<MPOptions>(static_cast<int>(f1) | static_cast<int>(f2));
}

constexpr MPOptions operator&(MPOptions f1, MPOptions f2) {
  return static_cast<MPOptions>(static_cast<int>(f1) & static_cast<int>(f2));
}

template <typename StreamType>
inline StreamType &operator<<(StreamType &s, MPOptions f) {
  s << "Used options: \n";
  if ((f & MPOptions::complex) != MPOptions::none)
    s << "	 - Complex spherical harmonics\n";
  if ((f & MPOptions::normalisation) != MPOptions::none)
    s << "	 - Normalisation is included\n";
  return s;
}
}  // namespace mutlipole_conv
#endif
