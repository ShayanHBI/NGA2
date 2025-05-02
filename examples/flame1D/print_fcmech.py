import os
import ARCANE.mechanisms as mechanisms
import ARCANE.custom.custom_kinetics as custom

mech = mechanisms.Mechanism("cti/gri30.cti")
custom.print_fortran(mech, "fcmech.f90", routine_name="fcmech", use="NGA")
os.system("mv fcmech.f90 src/fcmech.f90")