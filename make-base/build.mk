.PHONY: always

# A .mk file is built by running make on a target equal to that file's basename in the file's directory.
# The .mk file is then touched with a more recent modified date.
%.target.mk: always
	if ! make -C $$(dirname $@) -f $@ -q $@; then \
  		make -C $$(dirname $@) -f $@ $@ && touch $@; \
  	fi
