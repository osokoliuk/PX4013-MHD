cabal exec -- ghc --make \
    -threaded -dynamic -shared -fPIC -flink-rts \
    -g -O2 -rtsopts -isrc src/Lib.hs cbits/init.c \
    -package parallel -o libhsffi.so
