{-# LANGUAGE ExplicitNamespaces #-}
{-# LANGUAGE NamedFieldPuns #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RecordWildCards #-}

module Fourier where

{-
Module      : Fourier.hs
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
import Data.Complex (Complex (..), conjugate, imagPart, realPart)
import Data.List (transpose)

-- Some type definitions
type WindowFunction = Complex Double -> Complex Double
type SignalFunction = Complex Double -> Complex Double
type Index_i = Complex Double
type Index_j = Complex Double
type Index_k = Complex Double
type N = Complex Double
type T = Complex Double
type M = Complex Double
type T0 = Complex Double
type Ak = Complex Double
type OmegaT = Complex Double
type C2 = Complex Double
type C3 = Complex Double
type C4 = Complex Double

real :: Double -> Complex Double
real x = x :+ 0

average :: [Complex Double] -> Complex Double
average x = sum x / (real . fromIntegral . length $ x)

gaussian :: N -> T -> OmegaT -> Complex Double -> Complex Double
gaussian n t omegaT x = exp (-((x - n / 2) / (2 * t * omegaT)) ** 2)

windowFunction :: N -> T -> OmegaT -> Complex Double -> Complex Double
windowFunction n t omegaT x =
    gaussian t n omegaT x
        - ( gaussian t n omegaT (real (-1 / 2))
                * (gaussian t n omegaT (x + t) + gaussian t n omegaT (x - t))
          )
            / (gaussian t n omegaT (real (-1 / 2) + t) + gaussian t n omegaT (real (-1 / 2) - t))

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

calculateCumulants :: M -> N -> T -> OmegaT -> [Ak] -> C2
calculateCumulants m n t omegaT ak =
    let
        g = windowFunction n t omegaT

        c2 :: M -> [Ak] -> [Ak] -> C2
        c2 m x y = m / (m - 1) * (average (zipWith (\x y -> x * y) x y) - average x * average y)
     in
        n * (c2 m ak (conjugate <$> ak)) / (t * sum [g (j * t / n) * conjugate (g (j * t / n)) | j <- real <$> [0, 1 .. realPart n - 1]])

main :: IO ()
main = do
    let m = real 100
        n = real 100
        t = real 100
        omegaT = 0.14
        ak = [(\k -> fourierCoefficient t0 n t k (windowFunction n t omegaT) (cos)) <$> (real <$> [0, 1 .. 99]) | t0 <- (real <$> [0.0, realPart t .. realPart $ m * t])]
    print $ realPart <$> (\x -> calculateCumulants m n t omegaT x) <$> transpose ak
