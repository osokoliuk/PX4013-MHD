#!/usr/bin/env cabal
{- cabal:
build-depends:
  base,
  parallel
ghc-options: -main-is Fourier -O2 -threaded
-}
{-# LANGUAGE ExplicitNamespaces #-}
{-# LANGUAGE NamedFieldPuns #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RecordWildCards #-}

module Main where

{-
Module      : Main.hs
Description : Discrete Fourier methods module
Copyright   : (c) Oleksii Sokoliuk, 2026
License     : MIT
Maintainer  : oleksii.sokoliuk@mao.kiev.ua
Stability   : experimental
Portability : portable

A module that defines several useful functions computed in parallel,
namely Discrete Fourier Transform, Power Spectrum and higher order
polyspectra
-}

-- Module imports

import Control.Parallel.Strategies
import Data.Complex (Complex (..), conjugate, imagPart, magnitude, realPart)
import Data.List (elemIndex, minimumBy, transpose, zipWith4)
import qualified Data.Map as MP
import Data.Maybe (fromJust, fromMaybe)
import Data.Ord (comparing)
import Numeric.Tools.Integration
import System.Random

-- Some type definitions
type WindowFunction = Complex Double -> Complex Double
type SignalFunction = Complex Double -> Complex Double
type Index_k = Complex Double
type Index_l = Complex Double
type Index_p = Complex Double
type N = Complex Double
type T = Complex Double
type M = Complex Double
type T0 = Complex Double
type Ak = Complex Double
type Al = Complex Double
type Ap = Complex Double
type Akl = Complex Double
type Aklp = Complex Double
type OmegaT = Complex Double
type C2 = Complex Double
type C3 = Complex Double
type C4 = Complex Double
type S2 = Complex Double
type S3 = Complex Double
type S4 = Complex Double
type Time = Complex Double

-- Datatype definitions
data FourierCoefficientKind
    = DiscreteFourier
    | ContinuousFourier
    deriving (Eq, Show, Read)

-- Data type that describes the precision that you want to achieve with the
-- Gaussian quadrature integration, the only plausible choices are:
--  * 1e-6, 1e-7, 1e-8, 1e-9
--   (fast) ------> (slow)
data Precision
    = P6
    | P7
    | P8
    | P9
    deriving (Eq, Show, Read, Ord)

{- | Generate an integrator based on the given precision,
note that in this code we are only using Gaussian quadratures
as an integration method
-}
makeIntegrator :: Precision -> QuadParam
makeIntegrator precision =
    case precision of
        P6 -> QuadParam{quadPrecision = 1e-6, quadMaxIter = 20}
        P7 -> QuadParam{quadPrecision = 1e-7, quadMaxIter = 20}
        P8 -> QuadParam{quadPrecision = 1e-8, quadMaxIter = 20}
        P9 -> QuadParam{quadPrecision = 1e-9, quadMaxIter = 20}
        _ -> error "Incorrect precision given"

realIntegral :: QuadParam -> (Double -> Double) -> (Double, Double) -> Double
realIntegral integrator f (lo, hi) = fromMaybe 0 . quadRes $ quadRomberg integrator (lo, hi) f

complexIntegral :: QuadParam -> (Complex Double -> Complex Double) -> Complex Double -> Complex Double -> Complex Double
complexIntegral integrator f a b = r :+ i
  where
    r = realIntegral integrator realF (0, 1)
    i = realIntegral integrator imagF (0, 1)
    realF t = realPart (f (interpolate t)) -- or realF = realPart . f . interpolate
    imagF t = imagPart (f (interpolate t))
    interpolate t = a + (t :+ 0) * (b - a)

-- | Define a shortcut for a number with Im(x) = 0
real ::
    Double ->
    Complex Double
real x = x :+ 0

-- | Average value from an array of numbers
average ::
    [Complex Double] ->
    Complex Double
average x = sum x / (real . fromIntegral . length $ x)

twoAverage ::
    [Complex Double] ->
    [Complex Double] ->
    Complex Double
twoAverage x y = average (zipWith (\x y -> x * y) x y)

threeAverage ::
    [Complex Double] ->
    [Complex Double] ->
    [Complex Double] ->
    Complex Double
threeAverage x y z = average (zipWith3 (\x y z -> x * y * z) x y z)

fourAverage ::
    [Complex Double] ->
    [Complex Double] ->
    [Complex Double] ->
    [Complex Double] ->
    Complex Double
fourAverage x y z w = average (zipWith4 (\x y z w -> x * y * z * w) x y z w)

interpolate :: (Fractional a) => (a, a) -> (a, a) -> a -> a
interpolate (a, av) (b, bv) x = av + (x - a) * (bv - av) / (b - a)

mapLookup :: MP.Map Double Double -> Double -> Double
mapLookup m x =
    case (MP.lookupLE x m, MP.lookupGE x m) of
        (Just (a, av), Just (b, bv)) ->
            if a == b
                then av
                else interpolate (a, av) (b, bv) x
        (Nothing, Just (b, bv)) -> bv
        (Just (a, av), Nothing) -> av
        _ -> error "mapLookup"

closestValue :: (Ord a, Num a) => a -> [(a, b)] -> b
closestValue x m = snd $ minimumBy (comparing (abs . subtract x . fst)) m

-- | Gaussian function with an amplitude 1 and standard deviation t * sigmaT
gaussian ::
    N ->
    T ->
    OmegaT ->
    Complex Double ->
    Complex Double
gaussian n t sigmaT x = exp (-((x - n / 2) / (2 * t * sigmaT)) ** 2)

-- | Approximate confined Gaussian window function
windowFunction ::
    N ->
    T ->
    OmegaT ->
    Complex Double ->
    Complex Double
windowFunction n t sigmaT x =
    if realPart x < 0 || realPart x > realPart t
        then real 0.0
        else
            gaussian t n sigmaT x
                - ( gaussian t n sigmaT (real (-1 / 2))
                        * (gaussian t n sigmaT (x + t) + gaussian t n sigmaT (x - t))
                  )
                    / (gaussian t n sigmaT (real (-1 / 2) + t) + gaussian t n sigmaT (real (-1 / 2) - t))

-- | Check if a number can be represented by 2**x
isPowerOfTwo :: Int -> Bool
isPowerOfTwo 1 = True
isPowerOfTwo x
    | mod x 2 == 0 = isPowerOfTwo (div x 2)
    | otherwise = False

-- | Fast Fourier Transform (FFT) function, taken from https://github.com/sjappig/convhs
fftRaw :: (RealFloat a) => [a] -> Int -> Int -> [Complex a]
fftRaw _ 0 _ = []
fftRaw [] _ _ = []
fftRaw (x0 : _) 1 _ = [x0 :+ 0]
fftRaw x n s = zipWith (+) x1 x2 ++ (zipWith (-) x1 x2)
  where
    x1 = fftRaw x (div n 2) (2 * s)
    x2 = zipWith (*) [exp (0 :+ (-2 * pi * fromIntegral k / (fromIntegral n))) | k <- [0 .. ((div n 2) - 1)]] (fftRaw (drop s x) (div n 2) (2 * s))

-- | FFT helper function
fft :: (RealFloat a) => [a] -> [Complex a]
fft x
    | isPowerOfTwo n = fftRaw x n 1
    | otherwise = error "FFT works only for powers of two"
  where
    n = length x

-- | Inverse FFT function, taken from https://github.com/sjappig/convhs
ifftRaw :: (RealFloat a) => [Complex a] -> Int -> Int -> [Complex a]
ifftRaw _ 0 _ = []
ifftRaw [] _ _ = []
ifftRaw (x0 : _) 1 _ = [x0]
ifftRaw x n s = zipWith (+) x1 x2 ++ (zipWith (-) x1 x2)
  where
    x1 = ifftRaw x (div n 2) (2 * s)
    x2 = zipWith (*) [exp (0 :+ (2 * pi * fromIntegral k / (fromIntegral n))) | k <- [0 .. ((div n 2) - 1)]] (ifftRaw (drop s x) (div n 2) (2 * s))

-- | Inverse FFT helper function
ifft :: (RealFloat a) => [Complex a] -> [a]
ifft x
    | isPowerOfTwo n = [realPart v / (fromIntegral n) | v <- ifftRaw x n 1]
    | otherwise = error "IFFT works only for powers of two"
  where
    n = length x

-- | Function convolution via FFT, taken from https://github.com/sjappig/convhs
convFft :: (RealFloat a) => [a] -> [a] -> [a]
convFft h x
    | length x == length h = ifft (zipWith (*) (fft h) (fft x))
    | otherwise = error "Kernel and input data must have the same length"

{- | Calculate Fourier coefficients for a given function and a given window function
where the size of the window is T, containing N samples and starting from t = t0
with a wavenumber k
-}
fourierCoefficient ::
    FourierCoefficientKind ->
    Precision ->
    T0 ->
    N ->
    T ->
    [Index_k] ->
    WindowFunction ->
    SignalFunction ->
    [Ak]
fourierCoefficient kind prec t0 n t ks g f =
    let i = (0 :+ 1)
     in case kind of
            DiscreteFourier ->
                let
                    summand j k =
                        let x = j * t / n
                         in g x
                                * f (x - t0)
                                * exp (2 * pi * i * k * x / t)
                 in
                    (\x -> x * t / n)
                        <$> [ sum
                                [ summand j k
                                | j <- real <$> [0 .. realPart n - 1]
                                ]
                            | k <- ks
                            ]
            ContinuousFourier ->
                let ts = real <$> [realPart t0, realPart ((t - t0) / (n - 1) + t0) .. realPart (t0 + t)]
                    fs = realPart <$> (f <$> ts)
                    gs = realPart <$> (g <$> ts)
                 in fft (zipWith (\x y -> 2 * pi * x * y) fs gs)
            _ -> error "Incorrect Fourier coefficient kind specified"

{- | Calculate the values for the unbiased cumulants C2, C3, C4 and the
corresponding polyspectra using the formulas
given in 10.1016/j.dsp.2026.105893
-}
calculatePolyspectra ::
    M ->
    N ->
    T ->
    OmegaT ->
    [[Ak]] ->
    [Index_k] ->
    [Index_l] ->
    [Index_p] ->
    ([S2], [[S3]], [[S4]])
calculatePolyspectra m n t sigmaT ak k l p =
    let
        -- Create map for a fast lookup of Fourier coefficients
        akMap :: MP.Map Double [Ak]
        akMap = MP.fromList $ zip ([0, 1 .. realPart n - 1]) (transpose ak)

        -- Shortcut for a window function
        g = windowFunction n t sigmaT

        -- Second order cumulant
        c2 :: M -> [Ak] -> [Ak] -> C2
        c2 m x y =
            m / (m - 1) * (twoAverage x y - average x * average y)

        -- Third order cumulant
        c3 :: M -> [Ak] -> [Ak] -> [Ak] -> C3
        c3 m x y z =
            (m ** 2 / ((m - 1) * (m - 2)))
                * ( threeAverage x y z
                        - twoAverage x y * average z
                        - twoAverage x z * average y
                        - twoAverage y z * average x
                        + 2 * average x * average y * average z
                  )

        -- Fourth order cumulant
        c4 :: M -> [Ak] -> [Ak] -> [Ak] -> [Ak] -> C4
        c4 m x y z w =
            ( m ** 2
                / ((m - 1) * (m - 2) * (m - 3))
                * ((m + 1) * fourAverage x y z w)
                - (m + 1) * (threeAverage x y z * average w + threeAverage x y w * average z + threeAverage x z w * average y + threeAverage y z w * average x)
                - (m - 1) * (twoAverage x y * twoAverage z w + twoAverage x z * twoAverage y w + twoAverage x w * twoAverage y z)
                + 2 * m * (twoAverage x y * average z * average w + twoAverage x z * average y * average w + twoAverage x w * average y * average z + twoAverage y z * average x * average w + twoAverage y w * average x * average z + twoAverage z w * average x * average y)
                - 6 * m * average x * average y * average z * average w
            )

        -- Second order polyspectrum (i.e. power spectrum)
        s2 :: M -> N -> T -> [Ak] -> S2
        s2 m n t ak =
            (n * (c2 m ak (conjugate <$> ak)))
                / (t * sum [g (j * t / n) * conjugate (g (j * t / n)) | j <- real <$> [0, 1 .. realPart n - 1]])

        -- Third order polyspectrum (i.e. bispectrum)
        s3 :: M -> N -> T -> MP.Map Double [Ak] -> Index_k -> Index_l -> S3
        s3 m n t akMap k l =
            let ak = fromMaybe [real 0] $ MP.lookup (realPart k) akMap
                al = fromMaybe [real 0] $ MP.lookup (realPart l) akMap
                akl = fromMaybe [real 0] $ MP.lookup (realPart (k + l)) akMap
             in (n * (c3 m ak al (conjugate <$> akl)))
                    / (t * sum [g (j * t / n) ** 2 * conjugate (g (j * t / n)) | j <- real <$> [0, 1 .. realPart n - 1]])

        -- Fourth order polyspectrum (i.e. bispectrum)
        s4 :: M -> N -> T -> MP.Map Double [Ak] -> Index_k -> Index_l -> Index_p -> S4
        s4 m n t akMap k l p =
            let ak = fromMaybe [real 0] $ MP.lookup (realPart k) akMap
                al = fromMaybe [real 0] $ MP.lookup (realPart l) akMap
                ap = fromMaybe [real 0] $ MP.lookup (realPart p) akMap
                aklp = fromMaybe [real 0] $ MP.lookup (realPart (k + l + p)) akMap
             in (n * (c4 m ak al ap (conjugate <$> aklp)))
                    / (t * sum [g (j * t / n) ** 3 * conjugate (g (j * t / n)) | j <- real <$> [0, 1 .. realPart n - 1]])

        -- Slice of a trispectrum along the diagonal S4(w1,-w1,w2), retains most of the information
        -- while being two dimensional and much cheaper to calculate
        compactS4 :: M -> N -> T -> MP.Map Double [Ak] -> Index_k -> Index_l -> S4
        compactS4 m n t akMap k l =
            let ak = fromMaybe [real 0] $ MP.lookup (realPart k) akMap
                al = fromMaybe [real 0] $ MP.lookup (realPart l) akMap
             in (n * (c4 m ak (conjugate <$> ak) al (conjugate <$> al)))
                    / (t * sum [g (j * t / n) ** 2 * (conjugate (g (j * t / n))) ** 2 | j <- real <$> [0, 1 .. realPart n - 1]])

        -- Evaluate all of the polyspectra in parallel
        evaluatedS2 = (parMap rdeepseq (\x -> s2 m n t x) ak)
        evaluatedS3 = parMap rdeepseq (\y -> map (\x -> s3 m n t akMap x y) k) l
        evaluatedCompactS4 = parMap rdeepseq (\y -> map (\x -> compactS4 m n t akMap x y) k) l
     in
        (evaluatedS2, evaluatedS3, evaluatedCompactS4)

main :: IO ()
main = do
    let gaussF x = exp (-x ** 2)
        randomList :: Int -> [Double]
        randomList seed = randoms (mkStdGen seed) :: [Double]

        kind = DiscreteFourier
        prec = P6
        m = real 10
        n = real 4096
        t = n * 0.1
        sigmaT = 0.14
        k = real <$> [0, 1 .. realPart n - 1]
        l = k
        p = k
        tstart = 0

    contents <- readFile "noise.dat"

    let
        fs = read <$> words contents
        ts = [0.0, 0.1 .. (0.1 * fromIntegral (length fs))]

        f x = mapLookup (MP.fromList $ zip ts fs) (realPart x)
        g x = windowFunction n t sigmaT x

        ak = parMap rpar (\t0 -> fourierCoefficient kind prec t0 n t k g (real . f)) (real <$> [tstart, tstart + realPart t .. realPart $ m * t])
        ak' = parMap rdeepseq (\t0 -> fourierCoefficient kind prec t0 n t k g (real . f)) (real <$> [tstart + realPart t / 2, tstart + realPart (3 * t / 2) .. realPart $ m * t])
        (s2', s3', s4') = (calculatePolyspectra m n t sigmaT (transpose ak) k l p)
        (s2'', s3'', s4'') = (calculatePolyspectra m n t sigmaT (transpose ak') k l p)
        s2 = zipWith (\x y -> (x + y) / 2) s2' s2''

    -- print $ (fmap . fmap) magnitude s4
    print $ realPart <$> s2
    -- print $ (fmap . fmap) magnitude s3
    print $ (fmap . fmap) realPart ak
    print ""
