all_drivers :
	cd drivers ; make all_drivers

lib :
	cd srcMPI  ; make BridgeMPI.a
	cd srcMT   ; make BridgeMT.a
	cd srcST   ; make Bridge.a

clean :
	cd srcMPI  ; make clean
	cd srcMT   ; make clean
	cd srcST   ; make clean
	cd drivers ; make clean
	cd doc     ; make clean
