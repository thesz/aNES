-- |ANES.hs
--
-- Approximate Natural Evolution Strategies.
--
-- Natural Evolution Strategy at Wikipedia: https://en.wikipedia.org/wiki/Natural_evolution_strategy
--
-- The code here follows spirit of NES in the algorithm and uses
-- approximation of square root of covariance matrix with
-- sum of diagonal and outer product matrices.
--
-- Tests indicate it kinda sorta works.
--
-- I am not sure about Fisher Information Matrix multiplication, I'll sort that
-- later. And I should have gradient of distribution in closed form, current
-- implementation is way too slow.
--
-- Copyright (C) 2020 Serguey Zefirov.

{-# OPTIONS -Wno-tabs #-}

{-# LANGUAGE FlexibleContexts, GeneralizedNewtypeDeriving
           , MultiParamTypeClasses, RecordWildCards
           , StandaloneDeriving #-}

module ANES where

import Control.Monad
import Control.Monad.State

import qualified Data.List as List

import System.Random

import qualified Numeric.LinearAlgebra as LA

import Debug.Trace

checkAsserts :: Bool
checkAsserts = True

type FLOAT = Double

type Vec = LA.Vector FLOAT

vSum :: Vec -> FLOAT
vSum = LA.sumElements

vProduct :: Vec -> FLOAT
vProduct = LA.prodElements

vMap :: (FLOAT -> FLOAT) -> Vec -> Vec
vMap f a = LA.cmap f a

vInner :: Vec -> Vec -> FLOAT
vInner a b = LA.sumElements $ a * b

vScale :: FLOAT -> Vec -> Vec
vScale a b = LA.cmap (* a) b

vL2, vL1 :: Vec -> FLOAT
vL2 a = vInner a a
vL1 a = vSum $ LA.cmap abs a

vFromList :: [FLOAT] -> Vec
vFromList = LA.fromList

vToList :: Vec -> [FLOAT]
vToList a = LA.toList a

vLength :: Vec -> Int
vLength a = LA.size a

-- |A=diag(D) + UV^T.
--
-- A is equal to diagonal matrix D plus outer product
-- of V and U.
data A = A { aD, aU, aV :: Vec } deriving (Show)

invertA :: A -> A
invertA (A d u v) = A d' u' v'
	where
		d' = vMap (1/) d
		g = vSum $ d' * u * v
		alpha = (-1) / (1 + g)
		v' = vScale alpha $ d' * v
		u' = d' * u

determinantA :: A -> FLOAT
determinantA (A d u v) = (1 + tr) * vProduct d
	where
		d' = vMap (1/) d
		tr = vSum $ d' * u * v

transposeA :: A -> A
transposeA (A d u v) = A d v u

matrMultA :: A -> Vec -> Vec
matrMultA (A d u v) w = dw + vScale vw u
	where
		dw = d * w
		vw = vInner v w

-- |Normal distribution parametrization.
--
-- S = AA^T, S^-1 = (A^-T)(A^-1).
-- PDF: P(x) = (2 * pi) ^ (- length mu / 2)
--           * (determinant(S)) ^ (-1/2)
--           * exp( (-1/2) * (x -mu ) (S^-1) (x - mu))
-- Juxtaposition above means matrix multiplication.
-- det(S) = det(A)^2
-- det(A) = (1 + trace(D^(-1)uv^T)) * det(D)
data N = N { nMu :: Vec, nA :: A } deriving (Show)

data DN = DN { dnDMu :: Vec, dnDA :: A } deriving (Show)

numDerivative :: (Vec -> FLOAT) -> Vec -> Vec
numDerivative f center = vScale (0.5/d) $ dps + dms
	where
		l = vLength center
		fc = f center
		d = 1e-4
		ds = map (vFromList . take l) $ take l $ iterate (0:) (d:repeat 0)
		dps = vFromList $ map ((\x -> x - fc) . f) $ map (center +) ds
		dms = vFromList $ map ((\x -> fc - x) . f) $ map (center -) ds

-- |Getting local gradient of a multivariate normal distribution
-- at point x.
--
-- Useful identities:
-- d/dS(det(S)) = (1/det(S)) (S^-1)^T
-- d/dA(S(A)^-1) = -(S(A)^-1)(d/dA(S(A)))(S(A)^-1)
-- d/dq(A(q)^TA(q)) = (dA/dq)A + A^TdA/dq
-- d/dD(D+uv^T) = I
-- d/du(D+uv^T) = v^T
-- d/dv(D+uv^T) = d/dv(u^Tv)^T = (d/dv(u^Tv))^T=u
dnAtX :: Vec -> N -> Vec
dnAtX x (N mu a) = numDerivative evalAtX v0
	where
		k = vLength x
		v0 = vFromList $ vToList mu ++ aToList a
		evalAtX v =	--m / sqrtDetS * exponent	-- makes for grad (    p(N,x))
				toExp - log (abs deta1) / 2	-- makes for grad (log(p(N,x))
			where
				m = (2 * pi) ** (negate $ fromIntegral k / 2)
				DN mu1 a1 = dnFromSplitted v
				a1' = invertA a1
				a1'T = transposeA a1'
				xmu = x - mu1
				a1'Txmu = matrMultA a1' xmu
				a1'Ta1'xmu = matrMultA a1'T a1'Txmu
				toExp = (/2) $ negate $ vInner xmu a1'Ta1'xmu
				exponent = exp toExp
				deta1 = determinantA a1
				sqrtDetS = abs deta1

dnFromSplitted :: Vec -> DN
dnFromSplitted v 
	| mod vl 4 /= 0 = error "invalid length at dnFromSplitted"
	| otherwise = DN mu a
	where
		vl = vLength v
		k = div vl 4
		mu = vFromList $ take k $ vToList v
		a = aFromList mu $ drop k $ vToList v

aToList :: A -> [FLOAT]
aToList (A d u v) = concatMap vToList [d,u,v]

aFromList :: Vec -> [FLOAT] -> A
aFromList v xs = A (vFromList dx) (vFromList ux) vx
	where
		l = vLength v
		(dx,xs') = splitAt l xs
		(ux,xs'') = splitAt l xs'
		vx = vFromList $ take l xs''

matrMultS :: A -> Vec -> Vec
matrMultS a = matrMultA (transposeA a) . matrMultA a

updateN :: FLOAT -> Vec -> N -> N
updateN nu dnVec n = N mu' a'
	where
		dn = dnFromSplitted dnVec
		mu' = nMu n + vScale nu (dnDMu dn)
		A d u v = nA n
		A dd du dv = dnDA dn
		a' = A
			(d + vScale nu dd)
			(u + vScale nu du)
			(v + vScale nu dv)

startN :: (Monad m, Draw m) => Vec -> m N
startN mu = do
	u <- drawStdNormalVec mu
	v <- drawStdNormalVec mu
	d <- drawStdNormalVec mu
	let	a = A d v u
	return $ N mu a

drawSample :: (MonadIO m, Draw m) => N -> m Vec
drawSample n = do
	x <- drawStdNormalVec $ nMu n
--	liftIO $ putStrLn $ "drawn stdnormal "++show x
	let	y = nMu n + matrMultA (nA n) x
--	liftIO $ putStrLn $ "transformed into "++show y
	return y

data Fisher = Fisher
	{ fisherN :: Double
	, gradients :: [Vec]
	}

conjugateGradient :: (Vec -> Vec) -> Vec -> Vec
conjugateGradient matMul b = cg r0 r0 x0
	where
		x0 = b
		r0 = b - matMul x0
		cg p r x
			| vL1 r < 1e-5 = x
			| vL1 r' < 1e-5 = x'
			| otherwise = cg p' r' x'
			where
				r2 = vInner r r
				ap = matMul p
				pAp = vInner p ap
				alpha = r2 / pAp
				x' = x + vScale alpha p
				r' = r - vScale alpha ap
				beta = vInner r' r' / r2
				p' = r' + vScale beta p

mulByFisherInverse :: Fisher -> Vec -> Vec
mulByFisherInverse m b = conjugateGradient multF b
	where
		fvecs = gradients m
		scale = 1 / fisherN m
		multF v = vScale (scale ^ 2) v + List.foldl1 (+)
			[vScale (scale * d) g | g <- fvecs, let d = vInner g v]

class Draw m where
	-- |Draw normally distributed (N(0, 1)) vector given another vector as a template.
	drawStdNormalVec :: Vec -> m Vec

-- |Approximate Natural Evolution Strategy.
--
-- Lambda is a population size. StartMu is an initial
-- averages. S is a state to evaluate against.
-- And eval is a State action to evaluate a sample. 
aNES :: (MonadIO m, Draw m) => Int -> FLOAT -> Vec -> (Vec -> m FLOAT) -> m (Maybe a) -> m a
aNES lambda nu startMu eval stop = do
	n <- startN startMu
	loop 0 n
	where
		ilambda = 1.0 / fromIntegral lambda
		loop i n = do
			liftIO $ putStrLn $ "iteration "++show i++", current n "++show n
			done <- stop
			case done of
				Just r -> return r
				_ -> contLoop i n
		drawAssess n = do
			sample <- drawSample n
			score <- eval sample
			let	grad = dnAtX sample n
			--liftIO $ putStrLn $ "sample "++show sample++", score "++show score++", grad "++show grad
			return (score, grad)
		contLoop i n = do
			scoredGrads <- fmap (zip (map ((^20) . (*ilambda)) [1..]) . map snd .
				List.sort) $ mapM (const $ drawAssess n) [1..lambda]
			let	matr = map snd scoredGrads
				grad = foldr1 (+) $ map (uncurry vScale) scoredGrads
				ngrad = mulByFisherInverse (Fisher (fromIntegral lambda) matr) grad
				n' = updateN nu ngrad n
			--liftIO $ putStrLn $ "matr "++show matr
			--liftIO $ putStrLn $ "grd "++show grad
			--liftIO $ putStrLn $ "ngrad "++show ngrad
			loop (i+1) n'

genNormals :: StdGen -> [FLOAT]
genNormals sg = transform rs
	where
		rs = map (\x -> x * 2 - 1) $ randoms sg :: [FLOAT]
		transform (u:v:uvs)
			| s == 0 || s >= 1 = transform uvs
			| otherwise = (ms * u) : (ms * v) : transform uvs
			where
				s = u^2 + v^2
				ms = sqrt $ negate $ 2 * log s / s


_test_cg1 = conjugateGradient mm (vFromList [1,-1.0])
	where
		mm ab = vFromList [b, a]
			where
				[a,b] = vToList ab

test_invertA1 = map (matrMultA a' . matrMultA a . vFromList) [[1,0,0], [0,1,0], [0,0,1],[1,1,1]]
	where
		a = A (vFromList [-1,2,3]) (vFromList [7,6,5]) (vFromList [11,13,17])
		a' = invertA a

test_numd1 = map (numDerivative eval . vFromList) [[1,1],[1,0],[0,1],[0,0]]
	where
		eval ab = (a*2 + b - 0.5)^2
			where
				[a,b] = vToList ab

data TEST1 = TEST1 { test1Best :: Vec, test1BestScore :: FLOAT, test1G :: [FLOAT] }
	deriving (Show)
newtype TEST1M a = TEST1M { runTEST1M :: StateT TEST1 IO a }
	deriving (Functor, Applicative, Monad, MonadIO)
deriving instance MonadState TEST1 TEST1M
instance Draw TEST1M where
	drawStdNormalVec v = do
		(ns, rest) <- fmap (splitAt (vLength v) . test1G) get
		modify $ \t1 -> t1 { test1G = rest }
		return $ vFromList ns
test_anes1 = fmap fst $ flip runStateT (TEST1 mu0 (-1e10) (genNormals $ mkStdGen 10)) $ runTEST1M $ do
	r <- aNES 10 0.5 (vFromList [10,10]) eval stop
	liftIO $ print r
	where
		mu0 = vFromList [10,10]
		eval v = do
			let score = 10 - (a*2 - 1)^2 - (b + 1)^2
			modify $ \t1@TEST1{..} -> if score > test1BestScore
				then TEST1 v score test1G
				else t1
			return score
			where
				[a,b] = vToList v
		stop :: TEST1M (Maybe Vec)
		stop = do
			s <- get
			let	bs = test1BestScore s
				bv = test1Best s
			--liftIO $ putStrLn $ "checking current best "++show bs++" for "++show bv
			--liftIO $ getLine
			if bs > 9.99999
				then do
					liftIO $ putStrLn $ "best score "++show bs++", specimen "++show bv
					return (Just bv)
				else return Nothing
t = test_anes1

