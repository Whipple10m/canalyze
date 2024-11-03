bary: bary.o bartim.o get_ephemeris.o earth.o unit.o loc.o sider.o tleap.o
	ar ruv libbary.a bary.o get_ephemeris.o bartim.o earth.o unit.o loc.o \
sider.o tleap.o
	ranlib libbary.a
#
bary.o: bary.c mystruc.h
	cc -c bary.c
#
get_ephemeris.o: get_ephemeris.c mystruc.h
	cc -c get_ephemeris.c
#
bartim.o: bartim.c mystruc.h
	cc -c bartim.c
#
earth.o: earth.c mystruc.h
	cc -c earth.c
#
unit.o: unit.c mystruc.h
	cc -c unit.c
#
loc.o: loc.c mystruc.h
	cc -c loc.c
#
sider.o: sider.c mystruc.h
	cc -c sider.c
#
tleap.o: tleap.c mystruc.h
	cc -c tleap.c
#
clean:
	rm -f *.o
#







