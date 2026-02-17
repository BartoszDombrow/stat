#' Szkielet: Bootstrap (Przedział Ufności i Estymatory)
#'
#' @description
#' Generuje kod do obliczenia estymatora, obciążenia i przedziału ufności metodą percentyli.
#' @export
kod_bootstrap <- function() {
  cat('
# --- BOOTSTRAP SZKIELET ---
# 1. DANE (Zmień to!)
dane <- ...
n <- length(dane)
theta_hat <- mean(dane) # <-- Zmień funkcję jeśli trzeba (median/var)

# 2. SYMULACJA
N <- 1000
boot <- numeric(N)

for(i in 1:N) {
  samp <- sample(dane, n, replace = TRUE)
  boot[i] <- mean(samp) # <-- Ta sama funkcja co wyżej!
}

# 3. WYNIKI
est_boot <- mean(boot)
bias <- mean(boot) - theta_hat
se <- sd(boot)

# Przedział ufności (Percentyle 95%)
alpha <- 0.05
k_low <- floor((N+1)*(alpha/2))
k_high <- floor((N+1)*(1-alpha/2))
sort_stat <- sort(boot)

cat("Estymator:", est_boot, "\\n")
cat("Bias:", bias, "\\n")
cat("CI 95%: [", sort_stat[k_low], ";", sort_stat[k_high], "]")
')
}

#' Szkielet: Jackknife (Obciążenie i Wariancja)
#'
#' @description
#' Generuje kod do metody Jackknife (pętla po n elementach).
#' @export
kod_jackknife <- function() {
  cat('
# --- JACKKNIFE SZKIELET ---
# 1. DANE
dane <- ...
n <- length(dane)
theta_hat <- mean(dane) # <-- Zmień funkcję

# 2. PĘTLA (tylko n razy!)
jack <- numeric(n)

for(i in 1:n) {
  samp <- dane[-i] # Usuwamy i-ty element
  jack[i] <- mean(samp) # <-- Zmień funkcję
}

# 3. WZORY JACKKNIFE
bias_jack <- (n - 1) * (mean(jack) - theta_hat)
var_jack <- ((n - 1) / n) * sum((jack - mean(jack))^2)
se_jack <- sqrt(var_jack)

cat("Bias:", bias_jack, "\\nSE:", se_jack)
')
}

#' Szkielet: Permutacje dla PAR (Zależne, np. przed/po)
#'
#' @description
#' Test permutacyjny dla dwóch pomiarów na tych samych obiektach.
#' Losowanie znaków +/-.
#' @export
kod_perm_pary <- function() {
  cat('
# --- PERMUTACJE PARY ---
x <- ... # Wektor 1
y <- ... # Wektor 2
diff <- x - y
sumo <- sum(diff)
n <- length(diff)

N <- 1000
csgo <- 0
stat <- numeric(n)

for(i in 1:N) {
  # Losowanie znaków
  stat <- ifelse(runif(n) < 0.5, diff, -diff)

  # WAŻNE: Sprawdź znak! (<= jeśli sumo ujemne, >= jeśli dodatnie)
  if(sum(stat) <= sumo) csgo <- csgo + 1
}

p_val <- 2 * csgo / N
cat("P-value:", p_val)
')
}

#' Szkielet: Permutacje NIEZALEŻNE (Dwie grupy)
#'
#' @description
#' Test dla dwóch różnych grup (wrzucanie do jednego worka).
#' @export
kod_perm_niezalezne <- function() {
  cat('
# --- PERMUTACJE NIEZALEŻNE ---
g1 <- ... # Grupa 1
g2 <- ... # Grupa 2 (testowana)
sumo <- sum(g2)
n <- length(g2)
kombi <- c(g1, g2) # Worek wspólny

N <- 1000
csgo <- 0

for(i in 1:N) {
  D <- sample(kombi, n) # Losujemy nową grupę
  # WAŻNE: Sprawdź znak! (<= jeśli sumo małe, >= jeśli duże)
  if(sum(D) <= sumo) csgo <- csgo + 1
}

p_val <- csgo / N
cat("P-value:", p_val)
')
}

#' Szkielet: Bootstrap PARAMETRYCZNY
#'
#' @description
#' Gdy znamy rozkład (Poisson, Wykładniczy, Normalny).
#' @export
kod_boot_param <- function() {
  cat('
# --- BOOTSTRAP PARAMETRYCZNY ---
x <- ...
n <- length(x)
# Parametr (np. lambda = mean(x) dla Poissona)
param <- mean(x)

set.seed(123)
B <- 5 # Ile próbek?
# Generowanie (rpois / rexp / rnorm)
boot_samples <- replicate(B, rpois(n, param))

# Statystyka dla każdej próbki (np. median)
boot_stats <- apply(boot_samples, 2, median)

bias <- mean(boot_stats) - median(x)
print(bias)
')
}
