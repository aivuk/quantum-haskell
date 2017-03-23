module Main where

import Numeric.GSL.Integration
import Complex
import System.IO
import Data.Fixed
import Data.Array
import Data.Maybe
import Control.Monad
import Data.Packed.Matrix
import Data.Packed.Vector
import Numeric.LinearAlgebra.Algorithms
import qualified Data.IntMap as M

data Ctree a b = End (a,b)
               | Branche (a,b) [Ctree a b]
    deriving (Show)

-- Given N boxes, return all possible combinations of P particles in these boxes.

c n p = g 1 0 
    where g i a | i == n     = [ End (n, p - a) ]
                | a == p - 1 = [ End (i, 1) ] ++ g (i + 1) a 
                | otherwise  = [ if r == p - a 
                                 then End (i,r) 
                                 else Branche (i,r) $ g (i+1) (a + r) | 
                                                         i <- [i..n-1],  
                                                         r <- [1..p-a] ] 
                                                            ++ [ End (n, p-a) ]

-- Convert combinations of Ctree do list format

cl [] = []
cl ((End p):l) = [p]:cl l
cl ((Branche p b):l) = (map (p:) $ cl b) ++ cl l

-- Create newstates aplying A+(i)A+(j)A(k)A(l) operator on state 'a'

newState i j k l n a = do cp i =<< cp j =<< dp k =<< dp l a 
    where
        cp i a = if a!i < n 
                 then Just $ a//[(i,a!i + 1)]
                 else Nothing

        dp i a = if a!i > 0
                 then Just $ a//[(i,a!i - 1)]
                 else Nothing

-- List all states that can have operator A+(i)A+(j)A(k)A(l) applied
 
possibleStates i j k l nb np = filter (isJust . newState i j k l np) . aStates nb

-- Autovalue of operator A+(i)A+(j)A(k)A(l) 

value i j k l n a = sqrt $ ((a!i + 1)*(a!j + 1)*a!k*a!l) / ((prod $ elems a')*(prod $ elems a))
    where a' = fromJust $ newState i j k l n a
          prod x = prod' x 1
          prod' [] a = a
          prod' (x:xs) a | x == 0 = prod' xs a
                         | otherwise = prod' xs a*x

sumValues i j k l n = sum . map (value i j k l n) 

-- Convert state from List to Array type

aStates nb = map (listArray (1,nb) (repeat 0) //) 

-- Total quasi-momentum of base state 'a'

qt a = qt' a 0 (min + 1) 
    where qt' a t i | i == max  = mod (t + (max - 1)*a!max) max
                    | otherwise = qt' a (t + (i-1)*a!i) (i + 1)
          (min,max) = bounds a

-- Group states by quasi-momentum

groupQt states = M.fromListWith (++) $ map (\x -> (qt x, [x])) states

-- Print grouped states

printGroups states = do mapM_ printG states
    where printG (n,s) = do putStrLn $ (show n) ++ ":"
                            mapM_ (putStrLn.show.elems) s
                            putStrLn ""

na = array (0, 4) [ (0, 1.81772),
                    (1, 1.75048),
                    (2, 0.998368),
                    (3, 0.111417),
                    (4, 1.19828) ]

a = array (0,4) [ (0, 1.94699),
                  (1, 1.98301),
                  (2, 2.29932),
                  (3, 2.49687),
                  (4, 2.22725) ]

n = 5 :: Int
nd = 5 :: Double

bloch phi w a b = exp(0 :+ a*phi)*exp(0 :+ (-wa*phy))*(cos(b*phy)*sin(b*pi/nd)*cos(wa*pi/nd) :+ sin(b*phy)*cos(b*pi/nd)*sin(wa*pi/nd))
    where phy = -pi/nd + mod' (phi + pi/nd) (2*pi/nd)
          wa = a - w

cbloch phi w a b = conjugate $ bloch phi w a b

toInt phi w i j k l = 1/((na!i)*(na!j)*(na!k)*(na!l))* 
                            cbloch phi w (fi i) (a!i) * cbloch phi w (fi j) (a!j) * 
                            bloch phi w (fi k) (a!k) * bloch phi w (fi l) (a!l)    
      where fi = fromIntegral 

-- Calculate matrix for each quasi-momentum

matFrom :: Int -> Int -> Array Int (Array Int Int) -> Double -> Array (Int,Int) Double
matFrom nb np gs c = matrix    
    where matrix = array ((1,1),(size,size)) [((i,j), f i j) | i <- [1..size], j <- [1..size]] 
          f i j | i > j =  matrix!(j,i)
                | otherwise = energy i j + c*(sum $ map (\x -> sq x (gs!j)*(im!x)) [ (oi,oj,ok,ol) | oi <- [1..nb], 
                                                                                                     oj <- [1..nb], 
                                                                                                     ok <- [1..nb], ol <- [1..nb], 
                                                                                                     newState oi oj ok ol np (gs!j) == Just (gs!i) ])
          sq (i,j,k,l) s = sqrt $ (fi $ (s!i) + 1)*(fi $ (s!j) + 1)*(fi $ s!k)*(fi $ s!l)
          energy i j | i == j = e (gs!j) 1 0 
                     | otherwise = 0
          e s i t | i == nb = t*(fi $ s!i)
                  | otherwise = e s (i+1) (t + (a!(i-1))*(fi $ s!i))
          size = snd $ bounds gs 
          im = intMatrix nb
          fi = fromIntegral

intMatrix n = array ((1,1,1,1),(n,n,n,n)) [ ((i,j,k,l), integrate (i,j,k,l) n) | i <- [1..n], j <- [1..n], k <- [1..n], l <- [1..n] ]

integrate (i,j,k,l) n | mod (i+j-k-l) n == 0 = fst $ integrateQAGS 1E-6 1000 (\x -> realPart $ toInt x 0.4 (i-1) (j-1) (k-1) (l-1)) 0 (2*pi)
                      | otherwise = 0

-- Create Matrix

genMatrix nb np k cst = matFrom nb np (listArray (1,length $ qmGroup M.! k) $ qmGroup M.! k) cst
    where qmGroup = groupQt $ aStates nb $ cl $ c nb np
          qmKeys = M.keys qmGroup    

-- Print eigenvalues and eigenstates to all states, group by quasi-momentum

printStates nb np cst file = do
    let qmGroup = groupQt $ aStates nb $ cl $ c nb np
        qmKeys = M.keys qmGroup    
 
        calEig n x = eig $ (n><n) $ elems x
        printEig (s,v) = do forM_ (zip (toList s) (toLists v)) $ \(e,v) -> do
                                 hPutStrLn file $ "Eigenvalue: " ++ (show e) ++ "\n"
                                 mapM_ (\n -> hPutStrLn file $ (show n) ++ " ") v
                                 hPutStrLn file "\n"

    forM_ qmKeys $ \k -> do
        let kSize = length $ qmGroup M.! k
        hPutStrLn file $ "quasi-momentum: " ++ (show k) ++ "\n"
        printEig $ calEig (kSize) $ matFrom nb np (listArray (1, kSize) $ qmGroup M.! k) cst

main = do
    file <- openFile "output" WriteMode
    printStates 5 10 1 file
    hClose file

