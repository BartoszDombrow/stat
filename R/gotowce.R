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

#' Szkielet: Permutacje dla PAR (Zależne - Wersja z dwiema pętlami)
#'
#' @description
#' Test permutacyjny dla par (np. klasyczne vs nowe).
#' Wersja zagnieżdżona (for w for), dokładnie jak na screenach z zajęć.
#' @export
kod_perm_pary <- function() {
  cat('
# --- PERMUTACJE PARY (Wersja z dwiema pętlami) ---
# 1. DANE
x <- ... # np. klasyczne
y <- ... # np. nowe
diff <- x - y
sumo <- sum(diff)
n <- length(diff)

# 2. PARAMETRY
N <- 1000
csgo <- 0
stat <- numeric(n)

# 3. SYMULACJA
for(i in 1:N) {

  # Pętla wewnętrzna (losowanie znaku dla każdego elementu z osobna)
  for(j in 1:n) {
    if(runif(1) < 0.5) {
      stat[j] <- diff[j]
    } else {
      stat[j] <- -diff[j]
    }
  }

  # Sprawdzanie warunku (ZNAK ZALEŻY OD SUMO!)
  # Jeśli sumo < 0, dajemy <=
  # Jeśli sumo > 0, dajemy >=
  if(sum(stat) <= sumo) {
    csgo <- csgo + 1
  }
}

# 4. WYNIK
p_val <- 2 * csgo / N
cat("Suma obserwowana (sumo):", sumo, "\\n")
cat("Licznik (csgo):", csgo, "\\n")
cat("P-value:", p_val, "\\n")
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

#' Szkielet: Bootstrap PARAMETRYCZNY (Poisson, Wykładniczy, Normalny)
#'
#' @description
#' Gdy znamy rozkład danych. Generuje kod z przykładami dla różnych rozkładów.
#' @export
kod_boot_param <- function() {
  cat('
# --- BOOTSTRAP PARAMETRYCZNY ---
x <- ... # Wstaw dane
n <- length(x)

# 1. ESTYMACJA PARAMETRÓW (Wybierz właściwy)
# --- Rozkład Poissona ---
param <- mean(x) # lambda_hat

# --- Rozkład Wykładniczy ---
# rate_hat <- 1/mean(x)

# --- Rozkład Normalny ---
# mu_hat <- mean(x)
# sd_hat <- sd(x)


# 2. GENEROWANIE PRÓBEK (Wybierz właściwy r...)
set.seed(123)
B <- 5 # Liczba próbek z polecenia

# Dla Poissona:
boot_samples <- replicate(B, rpois(n, param))

# Dla Wykładniczego:
# boot_samples <- replicate(B, rexp(n, rate_hat))

# Dla Normalnego:
# boot_samples <- replicate(B, rnorm(n, mu_hat, sd_hat))


# 3. STATYSTYKA I WYNIKI
# Liczymy statystykę dla każdej próbki (np. median, mean, var)
boot_stats <- apply(boot_samples, 2, median)

# Obciążenie (Bias) = (średnia z boot) - (statystyka z oryginału)
bias <- mean(boot_stats) - median(x)

cat("Statystyki z próbek:", boot_stats, "\\n")
cat("Obciążenie:", bias, "\\n")
')
}

#' Szkielet: Bootstrap dla KORELACJI (2 zmienne)
#'
#' @description
#' Generuje kod do zadania z macierzą kowariancji (rmvnorm) i korelacją.
#' UWAGA: Wymaga pakietu mvtnorm, jeśli generujesz dane.
#' @export
kod_boot_korelacja <- function() {
  cat('
# --- BOOTSTRAP KORELACJA (2D) ---
library(mvtnorm) # Odkomentuj jeśli generujesz dane!

# 1. PRZYGOTOWANIE DANYCH
# Opcja A: Generowanie (jak w zadaniu z rmvnorm)
# mu <- c(0,0)
# sigma <- diag(2)
# data <- rmvnorm(100, mean = mu, sigma = sigma)

# Opcja B: Masz gotowe wektory X i Y
# data <- data.frame(X = c(...), Y = c(...))

# --- Reszta jest automatyczna ---
n <- nrow(data)
index <- 1:n

# 2. PARAMETRY
N <- 999              # <-- WPISZ N (np. 999)
stat <- numeric(N)
conf_level <- 0.95    # <-- WPISZ POZIOM UFNOŚCI

# 3. PĘTLA (Symulacja)
for (i in 1:N) {
  # Klucz: Losujemy INDEKSY wierszy (żeby nie rozwalić par X-Y)
  idx_boot <- sample(index, n, replace = TRUE)

  # Tworzymy próbkę
  samp <- data[idx_boot, ]

  # Liczymy statystykę (tu: korelacja)
  stat[i] <- cor(samp[,1], samp[,2])
}

# 4. PRZEDZIAŁ UFNOŚCI (Percentyle)
sort_stat <- sort(stat)
alpha <- 1 - conf_level

# Automatyczne indeksy (dla N=999 to będzie 26 i 974)
k_low <- floor((N + 1) * (alpha / 2))
k_high <- floor((N + 1) * (1 - alpha / 2))

lower <- sort_stat[k_low]
upper <- sort_stat[k_high]
est <- mean(stat)

cat("Estymator:", est, "\\n")
cat("Przedział", conf_level*100, "%: [", lower, ";", upper, "]")
')
}

#' Szkielet: ANOVA - Blokowe Losowe (Kod ze zdjęcia)
#'
#' @description
#' Dokładne odwzorowanie kodu ze zdjęcia (Lista 6 Zadanie 1).
#' Liczy ręcznie sumy kwadratów (SST, SSA, SSB, SSE) i statystykę F.
#' @export
kod_anova_blok <- function() {
  cat('
# --- ANOVA: BLOKI LOSOWE (Kod ze zdjęcia) ---
# 1. WPROWADZANIE DANYCH
# Wpisz dane w c(...) i ustaw wymiary (4, 5) zgodnie z zadaniem
zad1 <- matrix(c(73, 68, 74, 71, 67,
                 73, 67, 75, 72, 70,
                 75, 68, 78, 73, 68,
                 73, 71, 75, 75, 69), 4, 5, byrow = T)

# 2. STAŁE
a <- nrow(zad1) # liczba wierszy (poziomów czynnika)
b <- ncol(zad1) # liczba kolumn (bloków)
N <- a * b      # liczebność całkowita

# Pomocnicza: "Suma całości do kwadratu / N" (Correction Factor)
CF <- (sum(zad1)^2)/N

# 3. SUMY KWADRATÓW (Wzory ze zdjęcia)
# SST - Całkowita zmienność
SST <- sum(zad1^2) - CF

# SSA - Wiersze (Czynnik badany)
# "suma kwadratów sum wierszowych dzielona przez liczbę kolumn"
SSA <- sum(rowSums(zad1)^2/b) - CF

# SSB - Bloki (Kolumny)
# "suma kwadratów sum kolumnowych dzielona przez liczbę wierszy"
SSB <- sum(colSums(zad1)^2/a) - CF

# SSE - Błąd (Reszta)
SSE <- SST - SSA - SSB

# 4. ŚREDNIE KWADRATY (MS)
MSA <- SSA / (a - 1)
MSB <- SSB / (b - 1)
MSE <- SSE / ((a - 1) * (b - 1))

# 5. TEST F
F_stat <- MSA / MSE

# 6. KWANTYL (Wartość krytyczna)
# qf(1-alfa, df1_czynnik, df2_blad)
F_kryt <- qf(0.95, (a - 1), (a - 1) * (b - 1))

# WYNIKI
cat("F-stat:", F_stat, "\\n")
cat("F-kryt (Kwantyl):", F_kryt, "\\n")

# Interpretacja
if(F_stat > F_kryt) {
  cat("Decyzja: ODRZUCAMY H0 (Istotne różnice)\\n")
} else {
  cat("Decyzja: BRAK PODSTAW do odrzucenia H0\\n")
}
')
}

#' Szkielet: ANOVA - Kwadrat Łaciński (Kod ze zdjęcia)
#'
#' @description
#' Dokładne odwzorowanie kodu ze zdjęcia (Zadanie 2 - Telewizory).
#' Uwzględnia macierz wartości (mat2) i macierz układu literek (literki).
#' @export
kod_anova_lacinski_orig <- function() {
  cat('
# --- ANOVA: KWADRAT ŁACIŃSKI (Oryginał) ---

# 1. DANE - WARTOŚCI (mat2)
# Wpisz dane kolumnami (tak jak na screenie: dół-góra, potem kolejna kolumna)
mat2 <- matrix(c(
  10, 7, 5, 10,   # Kolumna 1
  14, 18, 10, 10, # Kolumna 2
  7, 11, 11, 12,  # Kolumna 3
  8, 8, 9, 15     # Kolumna 4
), 4, 4)

# 2. UKŁAD LITEREK (literki)
# Zamień litery na liczby: A=1, B=2, C=3, D=4
# Wpisz układ kolumnami (tak jak wyżej)
literki <- matrix(c(
  4, 2, 1, 3,  # Kolumna 1 (D, B, A, C) -> (4, 2, 1, 3)
  3, 4, 2, 1,  # Kolumna 2
  1, 3, 4, 2,  # Kolumna 3
  2, 1, 3, 4   # Kolumna 4
), 4, 4)

# 3. OBLICZANIE SUM DLA LITEREK (Metoda ze zdjęcia "which")
SumA <- sum(mat2[which(literki==1)])
SumB <- sum(mat2[which(literki==2)])
SumC <- sum(mat2[which(literki==3)])
SumD <- sum(mat2[which(literki==4)])

# Wektor sum literek
suma_literek <- c(SumA, SumB, SumC, SumD)

# 4. STAŁE
N <- nrow(mat2) * ncol(mat2) # Liczebność całkowita (16)
a <- nrow(mat2)              # Wymiar kwadratu (4)

# 5. SUMY KWADRATÓW
# SST - Całkowita
SST <- sum(mat2^2) - (sum(mat2)^2/N)

# SSA - Litery/Kuracje (To co badamy - nazwa ze zdjęcia)
SSA <- sum(suma_literek^2/a) - (sum(mat2)^2/N)

# SSR - Wiersze (RowSums)
SSR <- sum(rowSums(mat2)^2/a) - (sum(mat2)^2/N)

# SSC - Kolumny (ColSums) - na screenie jako Operatorzy
SSC <- sum(colSums(mat2)^2/a) - (sum(mat2)^2/N)

# SSE - Błąd
SSE <- SST - SSA - SSR - SSC

# 6. ŚREDNIE KWADRATY I TEST F
MSA <- SSA / (a - 1)
# Stopnie swobody błędu dla kwadratu łacińskiego: (a-2)(a-1)
MSE <- SSE / ((a - 2) * (a - 1))

F0 <- MSA / MSE

# Wartość krytyczna (Kwantyl)
F_kryt <- qf(0.95, (a - 1), (a - 2) * (a - 1))

# WYNIKI
cat("Suma literek:", suma_literek, "\\n")
cat("SST:", SST, " SSA:", SSA, " SSE:", SSE, "\\n")
cat("F-stat:", F0, "\\n")
cat("F-kryt:", F_kryt, "\\n")

if(F0 > F_kryt) {
  cat("Decyzja: ODRZUCAMY H0 (Metody różnią się istotnie)\\n")
} else {
  cat("Decyzja: BRAK PODSTAW do odrzucenia H0\\n")
}
')
}


#' Szkielet: ANOVA - Kwadrat Grecko-Łaciński (Kod ze zdjęcia)
#'
#' @description
#' Dokładne odwzorowanie kodu ze zdjęcia (Zadanie 3 - Stanowisko pracy).
#' Obsługuje cztery źródła zmienności: Wiersze, Kolumny, Litery łacińskie i Litery greckie.
#' @export
kod_anova_grecko_lacinski <- function() {
  cat('
# --- ANOVA: KWADRAT GRECKO-ŁACIŃSKI ---

# 1. DANE - WARTOŚCI (mat3)
# Wpisz dane kolumnami (z góry na dół, potem kolejna kolumna)
mat3 <- matrix(c(
  11, 8, 9, 9,    # Kolumna 1
  10, 12, 11, 8,  # Kolumna 2
  14, 10, 7, 18,  # Kolumna 3
  8, 12, 15, 6    # Kolumna 4
), 4, 4)

# 2. UKŁAD LITER ŁACIŃSKICH (A=1, B=2, C=3, D=4)
literki_latino <- matrix(c(
  3, 1, 2, 4,     # Kolumna 1
  1, 3, 4, 2,     # Kolumna 2
  4, 2, 1, 3,     # Kolumna 3
  2, 4, 3, 1      # Kolumna 4
), 4, 4)

# 3. UKŁAD LITER GRECKICH (alpha=1, beta=2, gamma=3, delta=4)
literki_grekoraki <- matrix(c(
  2, 1, 3, 4,     # Kolumna 1
  3, 4, 2, 1,     # Kolumna 2
  4, 3, 1, 2,     # Kolumna 3
  1, 2, 4, 3      # Kolumna 4
), 4, 4)

# 4. SUMY DLA POSZCZEGÓLNYCH LITER (Metoda which)
# Sumy dla czynnika łacińskiego (np. Metody montażu)
latino_sumy <- c(
  sum(mat3[which(literki_latino==1)]),
  sum(mat3[which(literki_latino==2)]),
  sum(mat3[which(literki_latino==3)]),
  sum(mat3[which(literki_latino==4)])
)

# Sumy dla czynnika greckiego (np. Stanowiska pracy)
greko_sumy <- c(
  sum(mat3[which(literki_grekoraki==1)]),
  sum(mat3[which(literki_grekoraki==2)]),
  sum(mat3[which(literki_grekoraki==3)]),
  sum(mat3[which(literki_grekoraki==4)])
)

# 5. STAŁE
p <- nrow(mat3)         # Wymiar kwadratu (4)
N <- nrow(mat3) * ncol(mat3) # 16

# 6. SUMY KWADRATÓW (SS)
SST <- sum(mat3^2) - (sum(mat3)^2 / N)

# SSL - Łacińskie (np. Metody)
SSL <- sum(latino_sumy^2 / p) - (sum(mat3)^2 / N)

# SSG - Greckie (np. Stanowiska)
SSG <- sum(greko_sumy^2 / p) - (sum(mat3)^2 / N)

# SSR - Wiersze
SSR <- sum(rowSums(mat3)^2 / p) - (sum(mat3)^2 / N)

# SSC - Kolumny
SSC <- sum(colSums(mat3)^2 / p) - (sum(mat3)^2 / N)

# SSE - Błąd
SSE <- SST - SSL - SSG - SSR - SSC

# 7. ŚREDNIE KWADRATY I TEST F
# Testujemy czynnik grecki (tak jak na screenie: MSG / MSE)
MSG <- SSG / (p - 1)
MSE <- SSE / ((p - 3) * (p - 1)) # Stopnie swobody: (p-3)(p-1)

F0 <- MSG / MSE

# Wartość krytyczna
F_kryt <- qf(0.95, (p - 1), (p - 3) * (p - 1))

# WYNIKI
cat("SST:", SST, " SSL:", SSL, " SSG:", SSG, " SSE:", SSE, "\\n")
cat("F-stat (dla czynnika greckiego):", F0, "\\n")
cat("F-kryt:", F_kryt, "\\n")

if(F0 > F_kryt) {
  cat("Decyzja: ODRZUCAMY H0 (Czynnik ma istotny wpływ)\\n")
} else {
  cat("Decyzja: PRZYJMUJEMY H0 (Brak istotnego wpływu)\\n")
}
')
}

#' Szkielet: Bootstrap RESZTOWY (Wersja manualna z wektorami e1, e2, e3)
#' @export
kod_boot_regresja_manual <- function() {
  cat('
# --- BOOTSTRAP RESZTOWY (WERSJA MANUALNA) ---
# 1. DANE I MODEL BAZOWY
X <- c(4.7, 3.2, 5.3, 4.9, 2.9)
Y <- c(9.5, 2.5, 15.2, 9.2, 2.9)
eps <- c(-0.42, -0.52, 2.52, -1.64, 1.26)

# Parametry podane w zadaniu (y = B0 + B1*x)
B0 <- -11.7
B1 <- 4.6

# 2. WYZNACZENIE WARTOSCI TEORETYCZNYCH
Y_est <- B0 + (B1 * X)

# 3. LOSOWANIE WEKTORÓW RESZT (Zgodnie z indeksami w poleceniu)
e1 <- c(eps[2], eps[2], eps[5], eps[3], eps[5])
e2 <- c(eps[1], eps[2], eps[3], eps[5], eps[3])
e3 <- c(eps[4], eps[2], eps[3], eps[5], eps[4])

# 4. NOWE Y I MODELE
Y1 <- Y_est + e1
coef1 <- lm(Y1 ~ X)$coef

Y2 <- Y_est + e2
coef2 <- lm(Y2 ~ X)$coef

Y3 <- Y_est + e3
coef3 <- lm(Y3 ~ X)$coef

# 5. WYNIK KOŃCOWY (ŚREDNIA)
wyniki <- rbind(coef1, coef2, coef3)
print(wyniki)

cat("\\nŚrednie parametry (B0 i B1):\\n")
colMeans(wyniki)
')
}
