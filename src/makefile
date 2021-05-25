CC=gcc
CFLAGS=-Wall -std=gnu99

TARGET=GoJumpLca
SRCS =dfs.c father.c init.c LCA.c save_plot.c
OBJS = $(SRCS:.c=.o)

all:$(TARGET)

$(TARGET):libgolca.so
	$(CC) -o $@ main.c -L. -lgolca -I../include

libgolca.so:$(OBJS)
	$(CC) -shared -o $@ $^

%.o:%.c
	$(CC) $(CFLAGS) -fPIC -c -o $@ $<

.PHONY: clean move
clean:
	rm -rf $(TARGET) $(OBJS)

move:
	mv libgolca.so ../lib
	mv GoJumpLca ../bin
	mv *.o ../shared
