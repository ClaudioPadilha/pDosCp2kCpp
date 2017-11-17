CC = g++
OFLAG = -O2
CFLAGS = -std=c++11 -fpermissive $(OFLAG)
OBJ = dos_cp2k_v_1_1_0.cpp
OUT = dos_cp2k.x
T = test.cpp
O = test.x

obj: $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) -o $(OUT)
clean:
	rm -rf $(OUT)
test: $(T)
	$(CC) $(CFLAGS) $(T) -o $(O)
