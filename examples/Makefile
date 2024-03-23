.SUFFIXES:
.SUFFIXES: .pyt -calibrate.py -apply.py

.pyt-apply.py:
	tmpfile=`mktemp pyt2py_tempXXXXXX.py`; \
	./pyt2py.py "$<" > "$${tmpfile}"; \
	python3 "$${tmpfile}"; \
	rm -f "$${tmpfile}"

all: 1x1-corrected.s1p 1x1m-corrected.s1p 2x1-corrected.s2p \
	2x2ab-corrected.s2p TRL-corrected.s2p

1x1-corrected.s1p: 1x1-apply.py
	python 1x1-calibrate.py
	python 1x1-apply.py
	python compare-same.py 1x1-expected.s1p 1x1-corrected.s1p

1x1m-corrected.s1p: 1x1m-apply.py
	python 1x1m-calibrate.py
	python 1x1m-apply.py
	python compare-same.py 1x1m-expected.s1p 1x1m-corrected.s1p

2x1-corrected.s2p: 2x1-apply.py
	python 2x1-calibrate.py
	python 2x1-apply.py
	python compare-same.py 2x1-expected.s2p 2x1-corrected.s2p

2x2ab-corrected.s2p: 2x2ab-apply.py
	python 2x2ab-calibrate.py
	python 2x2ab-apply.py
	python compare-same.py 2x2ab-expected.s2p 2x2ab-corrected.s2p

TRL-corrected.s2p: TRL-apply.py
	python TRL-calibrate.py
	python TRL-apply.py
	python compare-same.py TRL-expected.s2p TRL-corrected.s2p

clean: FRC
	rm -f *-calibrate.py *-apply.py *.vnacal \
		*-corrected.s[1234]p *-expected.s[1234]p *-measured.s[1234]p \
		short.s1p open.s1p load.s1p

FRC:
