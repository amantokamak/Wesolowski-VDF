package main

import (
	"crypto/sha256"
	"fmt"
	"math/big"
	"os"
	"strconv"
	"sync"
	"time"
)

// Global variables
var G *big.Int
var PrimeL int64
var T int64
var Kappa int64
var Gamma int64
var C []*big.Int

// Hash function HG
func HG(x *big.Int) *big.Int {
	hash := sha256.Sum256(x.Bytes())
	return new(big.Int).SetBytes(hash[:])
}

// Setup function to initialize parameters
func Setup(g, primeL, t, kappa, gamma int64) {
	G = big.NewInt(g)
	PrimeL = primeL
	T = t
	Kappa = kappa
	Gamma = gamma
	C = make([]*big.Int, (t/(kappa*gamma))+1)
	modulus := big.NewInt(PrimeL)
	// Precompute values for C
	for i := int64(0); i <= t/(kappa*gamma); i++ {
		exp := new(big.Int).Exp(big.NewInt(2), big.NewInt(kappa*gamma*i), nil)
		C[i] = new(big.Int).Exp(G, exp, modulus)
	}
}

// TrapdoorSK function with optimizations
func TrapdoorSK(x *big.Int, t int64, sk *big.Int) (*big.Int, *big.Int, time.Duration) {
	start := time.Now()

	g := HG(x)
	e := new(big.Int).Exp(big.NewInt(2), big.NewInt(t), sk)
	y := new(big.Int).Exp(g, e, nil)
	l := HPrime(g, y)
	r := new(big.Int).Mod(new(big.Int).Exp(big.NewInt(2), big.NewInt(t), nil), l)
	q := new(big.Int).ModInverse(l, sk)
	q.Mul(q, new(big.Int).Sub(e, r))
	q.Mod(q, sk)
	pi := new(big.Int).Exp(g, q, nil)

	duration := time.Since(start)
	return y, pi, duration
}

// HPrime hash function for l that ensures the output is a prime number
func HPrime(g, y *big.Int) *big.Int {
    hash := sha256.New()
    hash.Write(g.Bytes())
    hash.Write(y.Bytes())
    candidate := new(big.Int).SetBytes(hash.Sum(nil))

    // Ensure the candidate is a prime number
    iteration := 0
    if candidate.Bit(0) == 0 { // If candidate is even, make it odd
        candidate.Add(candidate, big.NewInt(1))
    }

    for !candidate.ProbablyPrime(20) {
        candidate.Add(candidate, big.NewInt(2)) // Increment by 2 to stay odd and find next prime
        iteration++
        if iteration > 100000 { // Add a safety check to prevent infinite loop
            fmt.Println("Failed to find a prime number after 100000 iterations")
            break
        }
    }
    fmt.Printf("Prime candidate after %d iterations: %s\n", iteration, candidate.String())
    return candidate
}

// Verify function
func Verify(x, y, pi *big.Int, t int64) bool {
	g := HG(x)
	l := HPrime(g, y)

	r := new(big.Int).Exp(big.NewInt(2), big.NewInt(t), nil)
	r.Mod(r, l)

	lhs := new(big.Int).Exp(pi, l, G)
	rhs := new(big.Int).Exp(g, r, G)
	rhs.Mul(rhs, y)
	rhs.Mod(rhs, G)

	return lhs.Cmp(rhs) == 0
}

// OptimizedEval function for evaluation with worker pool
func OptimizedEval(x *big.Int, t, kappa, gamma int64) (*big.Int, *big.Int, time.Duration) {
	start := time.Now()

	modulus := big.NewInt(PrimeL)
	y := new(big.Int).Set(HG(x)) // Use the hashed x as the base
	pi := new(big.Int).SetInt64(1)

	// Use channels to collect results from workers
	yChan := make(chan *big.Int, gamma)
	piChan := make(chan *big.Int, gamma)

	var wg sync.WaitGroup

	// Worker function for parallel processing
	worker := func(j int64) {
		defer wg.Done()
		localY := new(big.Int).SetInt64(1)
		localPi := new(big.Int).SetInt64(1)

		tExp := new(big.Int).Exp(big.NewInt(2), big.NewInt(t), nil)

		for i := int64(0); i <= t/(kappa*gamma); i++ {
			b := getBlock(tExp, i*gamma+j, kappa, PrimeL)
			if b != 0 {
				localY.Mul(localY, C[i])
				localY.Mod(localY, modulus)

				localPi.Mul(localPi, C[i])
				localPi.Mod(localPi, G)
			}
		}

		yChan <- localY
		piChan <- localPi
	}

	// Launch workers
	for j := gamma - 1; j >= 0; j-- {
		wg.Add(1)
		go worker(j)
	}

	// Aggregate results from workers
	go func() {
		wg.Wait()
		close(yChan)
		close(piChan)
	}()

	for localY := range yChan {
		y.Mul(y, localY)
		y.Mod(y, modulus)
	}

	for localPi := range piChan {
		pi.Mul(pi, localPi)
		pi.Mod(pi, G)
	}

	duration := time.Since(start)
	return y, pi, duration
}

// getBlock function to split and compute blocks
func getBlock(e *big.Int, i, kappa, primeL int64) int64 {
	block := new(big.Int)
	exp := new(big.Int).Exp(big.NewInt(2), big.NewInt(i*kappa), nil)
	mod := new(big.Int).Exp(big.NewInt(2), big.NewInt(kappa), nil)
	block.Mod(exp.Mul(exp, e), mod.Mul(mod, big.NewInt(primeL)))
	return block.Int64()
}

func main() {
	// Fetch command-line arguments
	if len(os.Args) < 3 {
		fmt.Println("Usage: go run main.go <x_value> <t_value>")
		return
	}

	xValueStr := os.Args[1]
	tValueStr := os.Args[2]

	x, success := new(big.Int).SetString(xValueStr, 10)
	if !success {
		fmt.Println("Invalid value for x:", xValueStr)
		return
	}

	t, err := strconv.ParseInt(tValueStr, 10, 64)
	if err != nil {
		fmt.Println("Invalid value for t:", tValueStr)
		return
	}

	// Example values for other parameters
	g := int64(2)
	primeL := int64(101)
	kappa := int64(2)
	gamma := int64(2)

	// Setup
	Setup(g, primeL, t, kappa, gamma)

	// Measure TrapdoorSK performance
	y, pi, durationTrapdoorSK := TrapdoorSK(x, t, big.NewInt(primeL))
	fmt.Printf("y: %s\nπ: %s\n", y.String(), pi.String())
	fmt.Printf("TrapdoorSK Time: %v\n", durationTrapdoorSK)

	// Verify y and π
	valid := Verify(x, y, pi, t)
	fmt.Printf("Verification: %v\n", valid)

	// Measure OptimizedEval performance
	yOpt, piOpt, durationOptimizedEval := OptimizedEval(x, t, kappa, gamma)
	fmt.Printf("Optimized y: %s\nOptimized π: %s\n", yOpt.String(), piOpt.String())
	fmt.Printf("OptimizedEval Time: %v\n", durationOptimizedEval)

	// Verify OptimizedEval result
	validOpt := Verify(x, yOpt, piOpt, t)
	fmt.Printf("OptimizedEval Verification: %v\n", validOpt)
}
