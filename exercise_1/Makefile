CC = gcc

CFLAGS = -lm
NAME = mpc_2d

$(NAME).exe:
	$(CC) $(NAME).c $(CFLAGS) -o $@

.PHONY: clean

clean:
	rm -f $(NAME).exe

.PHONY: run

run: $(NAME).exe
	./$<
