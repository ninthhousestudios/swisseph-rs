SWISSEPH_DIR ?= ../../swisseph
GOLDEN_DIR   = ../tests/golden-data
CC           = cc
CFLAGS       = -Wall -I$(SWISSEPH_DIR)
LDFLAGS      = -L$(SWISSEPH_DIR) -lswe -lm

golden_gen: golden_gen.c $(SWISSEPH_DIR)/libswe.a
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

bench_c: bench_c.c $(SWISSEPH_DIR)/libswe.a
	$(CC) $(CFLAGS) -O2 -o $@ $< $(LDFLAGS)

$(GOLDEN_DIR)/date.json: golden_gen
	mkdir -p $(GOLDEN_DIR)
	./golden_gen date > $@

$(GOLDEN_DIR)/math.json: golden_gen
	mkdir -p $(GOLDEN_DIR)
	./golden_gen math > $@

generate: $(GOLDEN_DIR)/date.json $(GOLDEN_DIR)/math.json

clean:
	rm -f golden_gen bench_c

.PHONY: generate clean
