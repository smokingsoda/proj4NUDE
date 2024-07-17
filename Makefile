CC=x86_64-linux-gnu-gcc
CFLAGS = -g -Wall -std=c99 -fopenmp -mavx -mfma -D_GNU_SOURCE
LDFLAGS = -fopenmp

# 可执行文件
TARGET=test_main

# 目标文件
OBJS=test_main.o dumbmatrix.o matrix.o

# 目标规则
$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS) $(LDFLAGS) 

# 目标文件规则
test_main.o: test_main.c dumbmatrix.h matrix.h
	$(CC) $(CFLAGS) -c test_main.c

dumbmatrix.o: dumbmatrix.c dumbmatrix.h
	$(CC) $(CFLAGS) -c dumbmatrix.c

matrix.o: matrix.c matrix.h
	$(CC) $(CFLAGS) -c matrix.c

# 清理规则
clean:
	rm -f $(TARGET)
	rm -f *.o
