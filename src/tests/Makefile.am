TESTS =
check_PROGRAMS =

TEST_EXTENSIONS = .py
PY_LOG_COMPILER = $(PYTHON_BIN)

if BUILD_SIMULATOR
if HAVE_PYTHON
TESTS += simulator-dists.py
endif
endif

# make sure automake includes the script in a distribution
dist_check_SCRIPTS = simulator-dists.py
