SWISSEPH_DIR ?= ../../swisseph
GOLDEN_DIR   = ../tests/golden-data
CC           = cc
CFLAGS       = -Wall -I$(SWISSEPH_DIR)
LDFLAGS      = -L$(SWISSEPH_DIR) -lswe -lm

golden_gen: golden_gen.c $(SWISSEPH_DIR)/libswe.a
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

$(GOLDEN_DIR)/date.json: golden_gen
	mkdir -p $(GOLDEN_DIR)
	./golden_gen date > $@

generate: $(GOLDEN_DIR)/date.json

clean:
	rm -f golden_gen

.PHONY: generate clean
