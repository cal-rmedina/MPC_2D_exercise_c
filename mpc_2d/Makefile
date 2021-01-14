CXX = gcc

NVFLAGS = -Wall -O3
NAME = mpc_2d

$(NAME).exe: $(NAME).c
	$(CXX) $(CFLAGS) $(NVFLAGS) $< -o $@

.PHONY: clean

clean:
	rm -f $(NAME).exe

.PHONY: run

run: $(NAME).exe
	./$<
