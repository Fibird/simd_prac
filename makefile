
# file dirs
SRC_DIR					= ./src
BUILD_DIR				= ./build
OBJ_DIR					= ./$(BUILD_DIR)/obj

# compilation option
CC						= gcc
CFLAGS					= -g -Wall -mavx2 -lpthread
LIBS					= -lpthread

SRCS					= $(wildcard $(SRC_DIR)/*.c)
OBJS					= $(SRCS:$(SRC_DIR)/%.c=$(OBJ_DIR)/%.o)
EXES					= $(OBJS:$(OBJ_DIR)/%.o=$(BUILD_DIR)/%)

# don't remove objs
.SECONDARY: $(OBJS)

all: $(EXES)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	@mkdir -p $(OBJ_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILD_DIR)/%: $(OBJ_DIR)/%.o
	@mkdir -p $(BUILD_DIR)
	$(CC) -o $@ $< $(LIBS)

.PHONY: all clean

clean:
	rm -rf $(BUILD_DIR)
