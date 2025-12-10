target = neoahrs
build_dir = build

cpp_sources = \
libs/neoahrs/neoahrs.cpp \
\
libs/magcal/magcal.cpp \
libs/magcal/magsac.cpp \
libs/magcal/magseg.cpp \
\
libs/k_vector/k_vector.cpp \
\
libs/madgwick/madgwick.cpp \
libs/madgwick/madmadgwick.cpp \
\
libs/butterworth/butterworth.cpp \
\
main.cpp

c_sources = \
libs/attitude/dcm.c \
libs/attitude/axan.c \
libs/attitude/quat.c \
libs/attitude/euler.c

includes = \
-I. \
-Ilibs/Eigen \
-Ilibs/magcal \
-Ilibs/neoahrs \
-Ilibs/madgwick \
-Ilibs/attitude \
-Ilibs/k_vector \
-Ilibs/rapidcsv \
-Ilibs/butterworth

defines =

libs =
lib_dir =

c_compiler = gcc
cpp_compiler = g++

linker_flags = $(lib_dir) $(libs) -g
compiler_flags = -g -Og -Wall -Wextra -lm $(includes) $(defines)

objects += $(addprefix $(build_dir)/,$(notdir $(c_sources:.c=.o)))
vpath %.c $(sort $(dir $(c_sources)))

objects += $(addprefix $(build_dir)/,$(notdir $(cpp_sources:.cpp=.o)))
vpath %.cpp $(sort $(dir $(cpp_sources)))

all: $(build_dir)/$(target)

$(build_dir)/$(target): $(objects) Makefile
	$(cpp_compiler) $(objects) $(linker_flags) -o $@

$(build_dir)/%.o: %.c Makefile | $(build_dir)
	$(c_compiler) -c $(compiler_flags) $< -o $@

$(build_dir)/%.o: %.cpp Makefile | $(build_dir)
	$(cpp_compiler) -c $(compiler_flags) $< -o $@

$(build_dir):
	@mkdir $(build_dir)

run: all
	@$(build_dir)/$(target)

see:
	octave plot/plot_output.m
