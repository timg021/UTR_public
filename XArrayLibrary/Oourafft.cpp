#include "OouraFft.h"

template<> xar::_eValueType xar::OouraFft<float>::GetValuetype() { return eXAFloat; }
template<> xar::_eValueType xar::OouraFft<double>::GetValuetype() { return eXADouble; }