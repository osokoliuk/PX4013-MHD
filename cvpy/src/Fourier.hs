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
import Data.Complex (Complex (..))

-- Some type definitions
type WindowFunction = Complex Double -> Complex Double
type SignalFunction = Complex Double -> Complex Double
type Index_i = Complex Double
type Index_j = Complex Double
type Index_k = Complex Double
type N = Complex Double
type T = Complex Double
type T0 = Complex Double
type Ak = Complex Double

fourierCoefficient ::
    (Enum N) =>
    T0 ->
    T ->
    N ->
    Index_k ->
    WindowFunction ->
    SignalFunction ->
    Ak
fourierCoefficient t0 t n k g f =
    let i = (0 :+ 1)
     in (t / n) * 
      sum [g (j * t / n) * f (j * t / n - t0) 
        * exp (2 * pi * i * k * j / n) 
          * exp (-2 * pi * i * k * t0 / t) | j <- [0, 1 .. n - 1]]

main = print "1"
