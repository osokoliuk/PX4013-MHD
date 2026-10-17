cabal exec -- ghc --make \
    -threaded -dynamic -shared -fPIC -flink-rts \
    -g -O0 -isrc src/Lib.hs cbits/init.c \
    -package parallel -o libhsffi.so
