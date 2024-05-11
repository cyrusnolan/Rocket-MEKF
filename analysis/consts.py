import numpy as np

freq = 10
h = 1 / freq

# massprop
m_dry = 25600.0
m_fuel = 10000.0
length = 20.0
radius = 1.85928
I1 = 0.25*(m_dry+m_fuel)*radius*radius+(1.0/12.0)*(m_dry+m_fuel)*(2*length*2*length);
I2 = 0.25*(m_dry+m_fuel)*radius*radius+(1.0/12.0)*(m_dry+m_fuel)*(2*length*2*length);
I3 = 0.5*(m_dry+m_fuel)*radius*radius

# environment
g = 9.807
g0 = 9.807

# propulsion
Isp = 311.0
T_min = 0.40 * 411000.0
T_max = 411000.0
