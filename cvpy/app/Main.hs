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
import Data.Complex (Complex (..), conjugate, imagPart, realPart)
import Data.List (transpose)
import qualified Data.Map as MP

-- Some type definitions
type WindowFunction = Complex Double -> Complex Double
type SignalFunction = Complex Double -> Complex Double
type Index_i = Complex Double
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
threeAverage x y z = average (zipWith3 (\x y z -> x * y * y) x y z)

-- | Gaussian function with an amplitude 1 and standard deviation t * omegaT
gaussian ::
    N ->
    T ->
    OmegaT ->
    Complex Double ->
    Complex Double
gaussian n t omegaT x = exp (-((x - n / 2) / (2 * t * omegaT)) ** 2)

-- | Approximate confined Gaussian window function
windowFunction ::
    N ->
    T ->
    OmegaT ->
    Complex Double ->
    Complex Double
windowFunction n t omegaT x =
    gaussian t n omegaT x
        - ( gaussian t n omegaT (real (-1 / 2))
                * (gaussian t n omegaT (x + t) + gaussian t n omegaT (x - t))
          )
            / (gaussian t n omegaT (real (-1 / 2) + t) + gaussian t n omegaT (real (-1 / 2) - t))

{- | Calculate Fourier coefficients for a given function and a given window function
where the size of the window is T, containing N samples and starting from t = t0
with a wavenumber k
-}
fourierCoefficient ::
    T0 ->
    N ->
    T ->
    Index_k ->
    WindowFunction ->
    SignalFunction ->
    Ak
fourierCoefficient t0 n t k g f =
    let i = (0 :+ 1)
     in (t / n)
            * sum
                [ g (j * t / n)
                    * f (j * t / n - t0)
                    * exp (2 * pi * i * k * j / n)
                    * exp (-2 * pi * i * k * t0 / t)
                | j <- real <$> [0, 1 .. realPart n - 1]
                ]

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
    [[Al]] ->
    [[Ap]] ->
    [[Akl]] ->
    [[Aklp]] ->
    [S2]
calculatePolyspectra m n t omegaT ak al ap akl aklp =
    let
        akMap :: MP.Map Index_k Ak
        akMap = MP.fromList $ zip [0, 1 .. realPart n - 1] ak

        g = windowFunction n t omegaT

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

        -- Second order polyspectrum (i.e. power spectrum)
        s2 :: M -> T -> N -> [Ak] -> S2
        s2 m t n ak =
            (n * (c2 m ak (conjugate <$> ak)))
                / (t * sum [g (j * t / n) * conjugate (g (j * t / n)) | j <- real <$> [0, 1 .. realPart n - 1]])

        -- Third order polyspectrum (i.e. bispectrum)
        s3 :: M -> T -> N -> MP.Map Index_k Ak -> [Index_k] -> [Index_l] -> S2
        s3 m t n akMap k l =
            let ak = fromMaybe 0 MP.lookup k akMap
                al = fromMaybe 0 MP.lookup l akMap
                akl = fromMaybe 0 MP.lookup (k + l) akMap
             in (n * (c3 m ak al (conjugate <$> akl)))
                    / (t * sum [g (j * t / n) ** 2 * conjugate (g (j * t / n)) | j <- real <$> [0, 1 .. realPart n - 1]])

        evaluatedS2 = (parMap rdeepseq (\x -> s2 m t n x) ak)
        evaluatedS3 = (parMap rdeepseq (\x -> s3 m t n x y))
     in
        evaluatedS2

main :: IO ()
main = do
    let m = real 100
        n = real 100
        t = real 100
        omegaT = 0.14
        ak = parMap rdeepseq (\t0 -> (\k -> fourierCoefficient t0 n t k (windowFunction n t omegaT) (cos)) <$> (real <$> [0, 1 .. 99])) (real <$> [0.0, realPart t .. realPart $ m * t])
    print $ realPart <$> (calculatePolyspectra m n t omegaT (transpose ak) (transpose ak) (transpose ak) (transpose ak) (transpose ak))
