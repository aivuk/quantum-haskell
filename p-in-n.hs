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

-- New version

dist_p_in_n p n | n == 1 = [[p]]
                | p == 0 = [replicate n 0]
                | otherwise = concatMap (\x -> map (x:) (dist_p_in_n (p - x) (n - 1))) [0..p]

d_p_in_n p p_min n | n == 1 = [[p]]
                   | p == 0 = [replicate n 0]
                   | otherwise = [x:l | x <- [p_min..p], l <- d_p_in_n (p - x) p_min (n - 1)]
