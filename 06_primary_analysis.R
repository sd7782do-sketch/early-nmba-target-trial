## =============================================================================
## 06_primary_analysis.R (analysis specification v10)
## Blocco neuromuscolare precoce nell'insufficienza respiratoria ipossiemica
## grave — analisi primaria clone-censor-weight con overlap weighting,
## imputazione multipla e bootstrap.  Derivazione: MIMIC-IV.
##
## Prerequisiti: OUT_DIR nell'ambiente; cohort_A_sev.parquet da 03 e 05.
## Se esiste cohort_A_sev_dt.parquet (con `deathtime` unita da admissions)
## viene usato automaticamente al suo posto.
##
##  [18] diag_mean: aggregate() scartava tutte le righe quando la troncatura e'
##       disattivata, perche' w_cens_cap e pct_rows_truncated sono interamente
##       NA e la formula usa na.omit.
##
## Modifiche rispetto alla v9
##  [17] deathtime incongruenti. Nel file cohort_A_sev_dt.parquet 13 stay hanno
##       deathtime precedente al time zero, fino a 20 ore prima: sono
##       incongruenze dei record, non del join (le date di deathtime e dod
##       coincidono in tutti i casi). Vengono riportati a valore mancante, cosi'
##       ricadono sul percorso di dod invece di essere ancorati al primo
##       intervallo. Il conteggio finisce nel manifest.
##
## Modifiche rispetto alla v8
##  [16] La troncatura al 99esimo percentile e' DISATTIVATA nel primario e
##       spostata fra le sensibilita'. Nel run v8 il cap nel braccio comparatore
##       cadeva a 1,614: la distribuzione di w_cens e' concentrata attorno a 1,
##       quindi quell'1% di righe sopra il cap conteneva l'intera riponderazione
##       per censura informativa. Troncandolo i pesi diventavano quasi uniformi,
##       il rischio del comparatore scendeva a 0,327 (sotto il valore senza
##       covariate orarie) e il bilanciamento sotto i pesi finali peggiorava da
##       0,080 a 0,171. Il criterio di scelta e' il bilanciamento, non l'ESS:
##       l'ESS alta con pesi quasi uniformi e' un sintomo, non un pregio.
##
## Modifiche rispetto alla v7
##  [15] Troncatura dei pesi di censura al quantile W_TRUNC_Q entro braccio.
##       Con le covariate orarie il modello di censura discrimina molto meglio,
##       e alcuni pazienti che sopravvivono senza iniziare pur avendo altissima
##       probabilita' oraria di iniziare prendono pesi estremi: nel run v7 il
##       massimo nel braccio comparatore saliva a 72-108 con 99esimo percentile
##       fermo a 4,5, e l'ESS crollava da 1.600 a 816. La troncatura taglia
##       quella coda. L'analisi senza troncatura resta fra le sensibilita'.
##
## Modifiche rispetto alla v6
##  [13] Il modello di censura copiava in g0 solo le colonne di PS_VARS, ma
##       plateau non e' in PS_VARS (e' la somma di driving pressure e PEEP e
##       viene escluso dal punteggio di propensione). g0[["plateau"]] era NULL
##       e ifelse() falliva. Ora vengono copiate anche le colonne di TV_USE.
##  [14] Valori mancanti nelle covariate tempo-varianti. Il plateau manca nel
##       23% dei pazienti al basale: glm avrebbe scartato quelle righe e
##       fitted() sarebbe tornato piu' corto di g0, disallineando i pesi senza
##       errore. Ora il ripiego e' a catena (valore orario, poi basale, poi
##       driving pressure + PEEP, che dopo imputazione e' sempre disponibile),
##       una variabile ancora incompleta viene esclusa dal modello di censura
##       invece di corrompere i pesi, e un controllo esplicito verifica che
##       fitted() abbia la lunghezza di g0.
##
## Modifiche rispetto alla v5
##  [11] Covariate tempo-varianti ridotte a P/F, PEEP e pressione di plateau.
##       Gli equivalenti noradrenalinici sono esclusi: nella griglia oraria il
##       5,2% dei valori supera 5 mcg/kg/min, con 99esimo percentile 23 e
##       massimo 749, cioe' lo stesso difetto di conversione delle unita' gia'
##       visto al basale, amplificato dalla somma dei segmenti sovrapposti. Il
##       filtro di range li sostituirebbe con valori riportati in avanti proprio
##       nei pazienti con supporto emodinamico piu' pesante. Restano nel modello
##       come covariata basale.
##  [12] Diagnostica delle covariate tempo-varianti: quota di valori orari
##       scartati dal filtro di range e quota di ore in cui il valore corrente
##       differisce dal basale, stampate e nel manifest.
##
## Modifiche rispetto alla v4
##  [9]  Decessi noti come sola data collocati alla FINE del giorno di morte
##       (24 h) invece che a meta' giornata. Con l'offset a 12 h i decessi
##       dello stesso giorno finivano prima del loro tempo reale rispetto a
##       time zero: 179 morti dentro il grace period contro le 44 della
##       convenzione a fine giornata, e 41 eventi con tempo negativo ancorati
##       al primo intervallo. Sensibilita' a 12 e 18 h, con il numero di
##       eventi ancorati riportato accanto a ciascuna.
##  [10] Due punteggi di propensione con ruoli distinti: quello del tilt
##       stimato sull'intera coorte analitica (conserva il bilanciamento
##       esatto dell'overlap weighting), quello dell'aderenza stimato nel risk
##       set alla scadenza del grace (denominatore corretto del peso di
##       censura). Nella v4 un unico punteggio ristretto al risk set pesava
##       una coorte diversa da quella su cui era stimato, e il bilanciamento
##       basale del lattato saliva a 0,109 con sovracorrezione del segno.
##
## Modifiche introdotte nella v4 e mantenute
##  [6] Unita' di analisi: restrizione al PRIMO episodio qualificante per
##      paziente, applicata PRIMA delle esclusioni. Il parquet contiene 6.370
##      stay per 6.112 soggetti; ricampionare righe non era ricampionare
##      pazienti, un decesso poteva entrare due volte, e un paziente
##      curarizzato in un passaggio precedente rientrava come nuovo
##      utilizzatore. La restrizione a monte evita di selezionare il secondo
##      episodio condizionatamente all'esclusione del primo.
##  [7] Descrittiva per finestra di inizio (0-6, 6-12, 12-24, >24 h, mai):
##      confronta la gravita' a T0 fra strati e la mortalita' osservata. Serve
##      a distinguere il confondimento basale da quello che nasce dopo T0.
##  [8] Pesi di censura dipendenti dal tempo, OPZIONALI: se esiste
##      initiation/hourly_grace.parquet con le covariate orarie della finestra
##      (stay_id, hour, pf, peep, plateau, ne_equiv), il modello di censura del
##      clone "no early" le usa al posto dei soli valori basali. La censura
##      artificiale e' determinata dall'inizio dell'infusione, che nasce spesso
##      da un deterioramento successivo a T0: con sole covariate basali quella
##      censura resta informativa. Senza il file lo script gira come la v3 e lo
##      dichiara nel manifest.
##
## Modifiche introdotte nella v3 e mantenute
##  [1] Tempo di morte: si usa `deathtime` (ora esatta) quando disponibile;
##      solo per i decessi noti come data si colloca l'evento a DOD_OFFSET_H
##      ore dall'inizio del giorno di morte (24 h nel primario), con 12 e 18 h
##      come sensibilita' esplicita sull'imputazione del tempo di morte.
##  [2] Il punteggio di propensione per l'aderenza del clone "early" e' stimato
##      nel risk set alla scadenza del grace (chi e' ancora vivo a quel punto),
##      non nell'intera coorte: chi muore durante il grace non ha avuto
##      occasione di iniziare e non appartiene al denominatore pertinente.
##  [3] Bilanciamento sotto i pesi FINALI clone-censor (risk set post-grace),
##      accanto a quello dei soli pesi di overlap basali. Secondo love plot.
##  [4] La stima a fine grace period NON e' un controllo di falsificazione:
##      il clone "no early" viene censurato all'inizio dell'infusione, quindi
##      le due strategie divergono gia' dentro la finestra e il contrasto non
##      e' nullo per costruzione. Etichette rinominate: resta come descrittiva.
##  [5] Didascalia della Figura S1 corretta: pooled binomial con link
##      complementary log-log, non "pooled logistic".
##
## Convenzioni mantenute dalle versioni precedenti
##  - peso del clone "early" = 1 nel grace, 1/PS dopo; clone "no early" pesato
##    con il prodotto cumulativo dell'azzardo di inizio, poi congelato
##  - il tilt di overlap e' identico sui due cloni; dopo il grace il peso del
##    clone early e' e(X)[1-e(X)]/p_adh(X), mentre quello del clone no early
##    usa la probabilita' cumulativa stimata dal modello orario di aderenza
##  - griglia oraria nel grace, poi bordi su giorni interi; orizzonte 28 giorni
##  - filtro di plausibilita' sulle covariate, valori fuori range imputati
##  - nessuna censura alla dimissione
## =============================================================================

suppressPackageStartupMessages({
    library(DBI); library(duckdb); library(splines)
    library(mice); library(ggplot2)
})
if (!exists("OUT_DIR")) stop("Definire OUT_DIR prima di eseguire lo script.")

## ---- configurazione ----------------------------------------------------------
QUICK       <- FALSE          # TRUE = M 3 / B 25, solo per verificare che giri
M_IMP       <- if (QUICK) 3  else 10
B_BOOT      <- if (QUICK) 25 else 200
GRACE_H     <- 12             # grace period primario, in ore
GRACE_SENS  <- c(6, 24)       # sensibilita' sul grace
HORIZON_D   <- 28             # orizzonte ESATTO, in giorni
LANDMARKS   <- c(7, 14, 28)
TRUNC_P     <- 0.02
DOD_OFFSET_H      <- 24       # [9] fine del giorno di morte: garantisce tempi
                              #     positivi e non richiede troncamenti
DOD_OFFSET_SENS   <- c(12, 18)#     sensibilita' sulla collocazione
FIRST_STAY_ONLY   <- TRUE     # [6] un solo episodio qualificante per paziente
TV_VARS           <- c("pf", "peep", "plateau")   # [8][11] tempo-varianti
W_TRUNC_Q         <- 0        # [15][16] troncatura dei pesi di censura; 0 = disattivata
W_TRUNC_SENS      <- 0.99     #          quantile usato nella sola sensibilita'
SEED        <- 20260915

set.seed(SEED)
out_dir <- file.path(OUT_DIR, "primary"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
fig_dir <- file.path(out_dir, "figures");  dir.create(fig_dir, showWarnings = FALSE)
chk_dir <- file.path(out_dir, "checks");   dir.create(chk_dir, showWarnings = FALSE)

## ---- registro degli avvisi ---------------------------------------------------
## L'avviso "successi non interi" (pesi non interi in famiglia binomiale) e'
## benigno: gli SE del glm non vengono usati, la varianza e' bootstrap.
## Tutto il resto, in particolare la mancata convergenza, viene conservato.
WARN_LOG <- new.env(); WARN_LOG$msgs <- character(0)
TV_STATE <- new.env(); TV_STATE$reported <- FALSE
quiet_glm <- function(expr) {
    withCallingHandlers(expr, warning = function(w) {
        msg <- conditionMessage(w)
        if (!grepl("non inter|non-integer", msg))
            WARN_LOG$msgs <- c(WARN_LOG$msgs, msg)
        invokeRestart("muffleWarning")
    })
}
ess <- function(w) sum(w)^2 / sum(w^2)

## =============================================================================
## 1. Coorte
## =============================================================================
con <- dbConnect(duckdb::duckdb(), dbdir = ":memory:")
dbExecute(con, "SET threads=6"); dbExecute(con, "SET preserve_insertion_order=false")

COH_DT <- file.path(OUT_DIR, "initiation", "cohort_A_sev_dt.parquet")
COH    <- if (file.exists(COH_DT)) COH_DT else
    file.path(OUT_DIR, "initiation", "cohort_A_sev.parquet")
if (!file.exists(COH)) stop("Manca ", COH, ": eseguire 03 e 05.")
cat("Coorte:", basename(COH), "\n")

raw <- dbGetQuery(con, sprintf("SELECT * FROM read_parquet('%s')",
                               gsub("\\\\", "/", normalizePath(COH))))

## [8] griglia oraria della finestra, se disponibile
TV_FILE <- file.path(OUT_DIR, "initiation", "hourly_grace.parquet")
TV_OK   <- file.exists(TV_FILE)
tv_raw  <- if (TV_OK) dbGetQuery(con, sprintf("SELECT * FROM read_parquet('%s')",
                                gsub("\\\\", "/", normalizePath(TV_FILE)))) else NULL
dbDisconnect(con, shutdown = TRUE)

## ---- [6] un solo episodio qualificante per paziente -------------------------
## Applicata PRIMA delle esclusioni: selezionare il secondo episodio
## condizionatamente all'esclusione del primo sarebbe una selezione informativa.
n_stay_all <- nrow(raw); n_subj_all <- length(unique(raw$subject_id))
if (FIRST_STAY_ONLY && "subject_id" %in% names(raw)) {
    raw <- raw[order(raw$subject_id, raw$time_zero), ]
    raw <- raw[!duplicated(raw$subject_id), ]
}
cat("Episodi qualificanti:", n_stay_all, "per", n_subj_all, "pazienti; trattenuti:",
    nrow(raw), "\n")

b0 <- raw[raw$prevalent_user == 0 & raw$sed_active_t0 > 0, ]

flow <- data.frame(
    step = c("Qualifying ICU stays",
             "Excluded: subsequent qualifying stay in the same patient",
             "One stay per patient",
             "Excluded: prevalent NMBA user at time zero",
             "New users",
             "Excluded: no sedative infusion at time zero",
             "Analytic cohort"),
    n = c(n_stay_all,
          n_stay_all - nrow(raw),
          nrow(raw),
          sum(raw$prevalent_user > 0),
          sum(raw$prevalent_user == 0),
          sum(raw$prevalent_user == 0 & raw$sed_active_t0 == 0),
          nrow(b0)))
print(flow, row.names = FALSE)
write.csv(flow, file.path(out_dir, "strobe_flow.csv"), row.names = FALSE)

b0$gender            <- factor(b0$gender)
b0$anchor_year_group <- factor(b0$anchor_year_group)
b0$h_init <- ifelse(is.na(b0$h_to_init), Inf, b0$h_to_init)

## ---- [1] tempo di morte ------------------------------------------------------
## deathtime da admissions ha risoluzione oraria e copre i decessi intraospedalieri.
## dod da patients ha risoluzione giornaliera ma copre anche i decessi dopo la
## dimissione: si usa solo per chi non ha deathtime, collocando l'evento a
## offset ore dall'inizio del giorno di morte.
HAS_DEATHTIME <- "deathtime" %in% names(b0)
if (!HAS_DEATHTIME)
    message("`deathtime` assente dal parquet: si ricade su dod per tutti i decessi. ",
            "Aggiungerla in 03/05 per la risoluzione oraria.")

set_death_time <- function(d, offset_h) {
    ht <- rep(NA_real_, nrow(d))
    if (HAS_DEATHTIME)
        ht <- as.numeric(difftime(as.POSIXct(d$deathtime),
                                  as.POSIXct(d$time_zero), units = "hours"))
    only_date <- is.na(ht) & !is.na(d$dod)
    if (any(only_date))
        ht[only_date] <- as.numeric(
            difftime(as.POSIXct(d$dod) + offset_h * 3600,
                     as.POSIXct(d$time_zero), units = "hours"))[only_date]
    ## un evento non puo' cadere a o prima di time zero: si ancora al primo
    ## intervallo invece di essere perso silenziosamente
    clipped <- which(!is.na(ht) & ht <= 0)
    if (length(clipped)) ht[clipped] <- 0.5
    attr(ht, "n_clipped")   <- length(clipped)
    attr(ht, "n_date_only") <- sum(only_date)
    d$h_to_death <- ht
    d
}
## [17] deathtime precedente o uguale a time zero: record incongruente
n_dt_bad <- 0L
if (HAS_DEATHTIME) {
    bad_dt <- which(!is.na(b0$deathtime) &
                    as.POSIXct(b0$deathtime) <= as.POSIXct(b0$time_zero))
    n_dt_bad <- length(bad_dt)
    if (n_dt_bad) b0$deathtime[bad_dt] <- NA
}
b0 <- set_death_time(b0, DOD_OFFSET_H)
n_date_only <- attr(b0$h_to_death, "n_date_only")
n_clipped   <- attr(b0$h_to_death, "n_clipped")
n_date_only_28 <- sum(
    (if (HAS_DEATHTIME) is.na(b0$deathtime) else TRUE) &
    !is.na(b0$dod) & !is.na(b0$h_to_death) &
    b0$h_to_death > 0 & b0$h_to_death <= HORIZON_D * 24
)
cat("\nDecessi a 28 g:", sum(!is.na(b0$h_to_death) & b0$h_to_death <= HORIZON_D*24),
    "| di cui noti come sola data:", n_date_only_28,
    "| eventi ancorati al primo intervallo:", n_clipped, "\n")

## dose di vasopressore: zero strutturale, non mancante
if (all(c("ne_equiv", "vaso_active") %in% names(b0)))
    b0$ne_equiv[b0$vaso_active == 0 & is.na(b0$ne_equiv)] <- 0

## volume corrente per kg di peso predetto, se l'altezza e' disponibile
add_vt_pbw <- function(d) {
    hc <- intersect(c("height", "height_cm", "height_first"), names(d))
    if (!length(hc) || !"vt" %in% names(d)) return(d)
    hgt <- d[[hc[1]]]; hgt[hgt < 120 | hgt > 220] <- NA
    pbw <- ifelse(d$gender == "F", 45.5, 50) + 0.91 * (hgt - 152.4)
    d$vt_pbw <- d$vt / pbw
    d
}
b0 <- add_vt_pbw(b0)
VT_VAR <- if ("vt_pbw" %in% names(b0)) "vt_pbw" else
    if ("vt" %in% names(b0)) "vt" else character(0)
if (identical(VT_VAR, "vt"))
    message("Altezza non disponibile: il volume corrente entra nel PS in mL.")

## ---- filtro di plausibilita' -------------------------------------------------
b0_raw <- b0                 # copia non filtrata, per la sensibilita'
PLAUS <- list(pf = c(20, 600), peep = c(0, 30), fio2 = c(21, 100), pao2 = c(20, 600),
              plateau = c(5, 60), driving_pressure = c(1, 45), vt = c(100, 1200),
              vt_pbw = c(1, 20), lactate = c(0.2, 30), creat = c(0.1, 20),
              plt = c(1, 1500), ne_equiv = c(0, 5), age = c(18, 100))
impl <- list()
for (v in intersect(names(PLAUS), names(b0))) {
    x <- b0[[v]]; lo <- PLAUS[[v]][1]; hi <- PLAUS[[v]][2]
    bad <- which(x < lo | x > hi)
    impl[[v]] <- data.frame(variable = v, n_set_to_na = length(bad),
                            pct_of_cohort = round(100 * length(bad) / nrow(b0), 2),
                            min_observed = round(min(x, na.rm = TRUE), 1),
                            max_observed = round(max(x, na.rm = TRUE), 1),
                            retained_range = paste(lo, hi, sep = " - "))
    if (length(bad)) b0[[v]][bad] <- NA
}
impl <- do.call(rbind, impl)
cat("\n=== Valori fuori dal range di analisi, posti a NA ===\n")
print(impl, row.names = FALSE)
write.csv(impl, file.path(chk_dir, "range_filtering.csv"), row.names = FALSE)

PS_VARS <- c("age", "gender", "anchor_year_group", "pf", "peep", "fio2",
             "vaso_active", "ne_equiv", "inotrope_active",
             "lactate", "creat", "plt", "driving_pressure", VT_VAR)
PS_VARS <- PS_VARS[PS_VARS %in% names(b0)]
PS_FORM <- as.formula(paste("treat ~", paste(PS_VARS, collapse = " + ")))

## ---- [8] matrici orarie delle covariate tempo-varianti ----------------------
## Una matrice per variabile, righe = stay_id della coorte, colonne = ora 1..24.
## Il valore basale riempie l'ora 1 e viene portato avanti (LOCF) dove l'ora
## successiva manca, cosi' il modello di censura non perde righe.
TV_USE <- character(0); TV_MAT <- list(); TV_ID <- b0$stay_id
if (TV_OK && "stay_id" %in% names(tv_raw) && "hour" %in% names(tv_raw)) {
    TV_USE <- intersect(TV_VARS, intersect(names(tv_raw), names(b0)))
    H <- 24
    for (v in TV_USE) {
        M <- matrix(NA_real_, nrow = length(TV_ID), ncol = H)
        i <- match(tv_raw$stay_id, TV_ID)
        j <- as.integer(tv_raw$hour)
        ok <- !is.na(i) & j >= 1 & j <= H
        M[cbind(i[ok], j[ok])] <- tv_raw[[v]][ok]
        M[, 1] <- ifelse(is.na(M[, 1]), b0[[v]], M[, 1])   # ancoraggio al basale
        for (k in 2:H) M[, k] <- ifelse(is.na(M[, k]), M[, k - 1], M[, k])
        if (v %in% names(PLAUS)) {                          # stesso filtro di range
            M[M < PLAUS[[v]][1] | M > PLAUS[[v]][2]] <- NA
            for (k in 2:H) M[, k] <- ifelse(is.na(M[, k]), M[, k - 1], M[, k])
        }
        TV_MAT[[v]] <- M                                    # NA residui: valore basale
    }
    cat("Covariate tempo-varianti nel modello di censura:",
        paste(TV_USE, collapse = ", "), "\n")
    ## [12] quanto le tempo-varianti si discostano davvero dal basale
    TV_DIAG <- do.call(rbind, lapply(TV_USE, function(v) {
        M <- TV_MAT[[v]][, seq_len(min(24, ceiling(GRACE_H))), drop = FALSE]
        base <- b0[[v]]
        data.frame(variable = v,
                   pct_hours_na = round(mean(is.na(M)), 4),
                   pct_hours_differs_from_baseline =
                       round(mean(!is.na(M) & !is.na(base) & abs(M - base) > 1e-8), 4))
    }))
    cat("\n=== Covariate tempo-varianti nella finestra ===\n")
    print(TV_DIAG, row.names = FALSE)
} else {
    TV_DIAG <- NULL
    message("hourly_grace.parquet assente: il modello di censura usa solo ",
            "covariate basali. La censura artificiale resta potenzialmente ",
            "informativa rispetto al deterioramento successivo a T0.")
}

## =============================================================================
## 2. Griglia person-time e clonazione
## =============================================================================
make_grid <- function(grace_h) {
    bnd <- unique(c(seq(0, grace_h, by = 1),
                    seq(ceiling(grace_h / 24) * 24, HORIZON_D * 24, by = 24)))
    if (tail(bnd, 1) < HORIZON_D * 24) bnd <- c(bnd, HORIZON_D * 24)
    data.frame(interval = seq_len(length(bnd) - 1),
               t_start  = head(bnd, -1),
               t_end    = tail(bnd, -1),
               in_grace = as.integer(tail(bnd, -1) <= grace_h))
}

## braccio 0 ("no early"): censurato NELL'INTERVALLO in cui l'infusione inizia,
##                         se l'inizio cade dentro il grace
## braccio 1 ("early"):    censurato alla fine del grace se non ha iniziato
## Un decesso conta per entrambi i cloni solo se nessuno dei due ha gia' deviato.
build_pt <- function(dat, grace_h, g) {
    n <- nrow(dat); K <- nrow(g)
    row <- rep(seq_len(n), each = K)
    gg  <- g[rep(seq_len(K), times = n), ]
    hd  <- dat$h_to_death[row]; hi <- dat$h_init[row]

    mk <- function(a) {
        ev <- as.integer(!is.na(hd) & hd > gg$t_start & hd <= gg$t_end)
        ce <- if (a == 0)
            as.integer(is.finite(hi) & hi <= grace_h & hi > gg$t_start & hi <= gg$t_end)
        else
            as.integer(gg$t_end == grace_h & hi > grace_h)
        ce[ev == 1] <- 0L
        data.frame(row = row, arm = a, interval = gg$interval,
                   in_grace = gg$in_grace,
                   width_d = (gg$t_end - gg$t_start) / 24,
                   t_mid_d = (gg$t_start + gg$t_end) / 48,
                   event = ev, cens = ce)
    }
    out <- rbind(mk(0), mk(1))
    out$key <- paste(out$row, out$arm, sep = "_")
    stop_i <- tapply(ifelse(out$event == 1 | out$cens == 1, out$interval, NA),
                     out$key, function(x) if (all(is.na(x))) Inf else min(x, na.rm = TRUE))
    out[out$interval <= stop_i[out$key], ]
}

## =============================================================================
## 3. Pesi
##    [2] Il modello di aderenza del clone "early" e' stimato fra i pazienti
##    ancora a rischio alla scadenza del grace period. Chi muore prima non ha
##    avuto occasione di iniziare l'infusione e non appartiene al denominatore
##    della probabilita' di aderenza; le predizioni sono poi estese a tutti.
## =============================================================================
build_weights <- function(dat, pt, grace_h, base_w, w_trunc = W_TRUNC_Q) {
    ## [10] due punteggi, due ruoli distinti.
    ## tilt: definisce la popolazione bersaglio. Stimato sull'INTERA coorte
    ##   analitica, cosi' l'overlap weighting conserva la proprieta' di
    ##   bilanciamento esatto, che vale solo se il modello e' stimato sullo
    ##   stesso campione a cui i pesi si applicano.
    ## aderenza: e' il denominatore del peso di censura del clone "early" ed e'
    ##   un parametro di disturbo. Stimato fra chi e' ancora a rischio alla
    ##   scadenza del grace: chi muore prima non ha avuto occasione di iniziare.
    ## Conseguenza da dichiarare: il peso finale del clone early e'
    ## e(X)(1-e(X))/e_adh(X), non esattamente (1-e(X)).
    ps_tilt <- fitted(quiet_glm(glm(PS_FORM, family = binomial, data = dat)))
    ps_tilt <- pmin(pmax(as.numeric(ps_tilt), TRUNC_P), 1 - TRUNC_P)

    at_risk <- is.na(dat$h_to_death) | dat$h_to_death > grace_h
    f_adh   <- quiet_glm(glm(PS_FORM, family = binomial, data = dat[at_risk, ]))
    ps_adh  <- pmin(pmax(as.numeric(predict(f_adh, newdata = dat, type = "response")),
                         TRUNC_P), 1 - TRUNC_P)

    pt$w_cens <- 1
    sel <- pt$arm == 0 & pt$in_grace == 1
    if (any(sel)) {
        g0 <- pt[sel, ]
        ## [13] anche le colonne tempo-varianti, che non sono tutte in PS_VARS
        for (v in unique(c(PS_VARS, TV_USE))) g0[[v]] <- dat[[v]][g0$row]
        ## [8] dove disponibile, il valore corrente sostituisce quello basale:
        ## la censura artificiale nasce dall'inizio dell'infusione, che segue
        ## spesso un deterioramento avvenuto dentro la finestra.
        cens_vars <- PS_VARS
        if (length(TV_USE)) {
            k  <- match(dat$stay_id[g0$row], TV_ID)
            hh <- pmin(pmax(as.integer(g0$interval), 1), 24)
            for (v in TV_USE) {
                base_v <- g0[[v]]
                ## [14] ripiego: plateau = driving pressure + PEEP, sempre
                ## disponibile dopo imputazione
                if (v == "plateau" && all(c("driving_pressure", "peep") %in% names(dat)))
                    base_v <- ifelse(is.na(base_v),
                                     dat$driving_pressure[g0$row] + dat$peep[g0$row], base_v)
                cur <- TV_MAT[[v]][cbind(k, hh)]
                g0[[paste0(v, "_t")]] <- ifelse(is.na(cur), base_v, cur)
            }
            keep <- TV_USE[vapply(TV_USE,
                                  function(v) !anyNA(g0[[paste0(v, "_t")]]), TRUE)]
            if (!TV_STATE$reported) {
                TV_STATE$reported <- TRUE
                if (length(keep) < length(TV_USE))
                    message("Escluse dal modello di censura per valori mancanti residui: ",
                            paste(setdiff(TV_USE, keep), collapse = ", "))
                cat("Modello di censura, covariate tempo-varianti effettive:",
                    if (length(keep)) paste(keep, collapse = ", ") else "nessuna", "\n")
            }
            cens_vars <- c(setdiff(PS_VARS, keep), paste0(keep, "_t"))
        }
        h0 <- if (sum(g0$cens) >= 10) {
            m0 <- quiet_glm(glm(as.formula(paste("cens ~ t_mid_d +",
                                                 paste(cens_vars, collapse = " + "))),
                                family = binomial, data = g0))
            ## [14] se glm avesse scartato righe, i pesi sarebbero disallineati
            if (length(fitted(m0)) != nrow(g0))
                stop("Modello di censura: fitted() ha lunghezza ", length(fitted(m0)),
                     " contro ", nrow(g0), " righe. Covariate con valori mancanti.")
            pmin(pmax(fitted(m0), 0), 1 - TRUNC_P)
        } else rep(0, nrow(g0))
        o    <- order(g0$key, g0$interval)
        cum  <- ave(1 - h0[o], g0$key[o], FUN = cumprod)
        g0$w <- NA_real_; g0$w[o] <- 1 / pmax(cum, TRUNC_P)
        pt$w_cens[sel] <- g0$w
        last0 <- tapply(g0$w, g0$key, function(x) x[length(x)])
        post  <- pt$arm == 0 & pt$in_grace == 0
        pt$w_cens[post] <- ifelse(is.na(last0[pt$key[post]]), 1, last0[pt$key[post]])
    }
    a1 <- pt$arm == 1
    pt$w_cens[a1] <- ifelse(pt$in_grace[a1] == 1, 1, 1 / ps_adh[pt$row[a1]])

    ## [15] troncatura della coda superiore del peso di censura, entro braccio.
    ## Si tronca w_cens e non il peso finale, perche' e' li' che nasce la coda:
    ## il tilt di overlap e' limitato per costruzione.
    trunc_info <- NULL
    if (!is.null(w_trunc) && w_trunc > 0 && w_trunc < 1) {
        rows <- list()
        for (a in 0:1) {
            sel_a <- which(pt$arm == a)
            cap   <- unname(quantile(pt$w_cens[sel_a], w_trunc, na.rm = TRUE))
            n_cut <- sum(pt$w_cens[sel_a] > cap, na.rm = TRUE)
            pt$w_cens[sel_a] <- pmin(pt$w_cens[sel_a], cap)
            rows[[length(rows) + 1]] <- data.frame(
                arm = a, w_cens_cap = round(cap, 3),
                pct_rows_truncated = round(n_cut / length(sel_a), 4))
        }
        trunc_info <- do.call(rbind, rows)
    }

    tilt <- switch(base_w, overlap = ps_tilt * (1 - ps_tilt), none = rep(1, length(ps_tilt)))
    pt$w <- pt$w_cens * tilt[pt$row]
    pt$w <- pt$w / mean(pt$w)      # normalizzazione globale: non altera i rapporti
    attr(pt, "ps")           <- ps_tilt
    attr(pt, "n_at_risk_ps") <- sum(at_risk)
    attr(pt, "trunc")        <- trunc_info
    pt
}

## =============================================================================
## 4. Stimatori
## =============================================================================
grid_time <- function(g) {
    tm <- data.frame(interval = g$interval,
                     t_mid_d  = (g$t_start + g$t_end) / 48,
                     width_d  = (g$t_end - g$t_start) / 24)
    tm$day_end <- cumsum(tm$width_d); tm
}

## azzardi discreti saturi con gli stessi pesi: versione non parametrica
sat_curve <- function(pt, g) {
    tm <- grid_time(g)
    r <- lapply(0:1, function(a) {
        s   <- pt[pt$arm == a, ]
        num <- tapply(s$w * s$event, s$interval, sum)
        den <- tapply(s$w, s$interval, sum)
        h <- rep(0, nrow(tm)); i <- as.integer(names(den))
        h[i] <- as.numeric(num / den)
        1 - cumprod(1 - h)
    })
    data.frame(day = tm$day_end, risk0 = r[[1]], risk1 = r[[2]], rd = r[[2]] - r[[1]])
}

estimate <- function(pt, g) {
    m <- quiet_glm(
        glm(event ~ arm + ns(t_mid_d, df = 3) + arm:ns(t_mid_d, df = 3) +
                offset(log(width_d)),
            family = binomial(link = "cloglog"), weights = w, data = pt))
    tm <- grid_time(g)
    cif <- lapply(c(0, 1), function(a) {
        nd <- data.frame(arm = a, t_mid_d = tm$t_mid_d, width_d = tm$width_d)
        h  <- unname(1 - exp(-exp(predict(m, newdata = nd, type = "link"))))
        1 - cumprod(1 - h)
    })
    r0 <- cif[[1]]; r1 <- cif[[2]]
    est <- c(risk0 = tail(r0, 1), risk1 = tail(r1, 1),
             rd = tail(r1, 1) - tail(r0, 1), rr = tail(r1, 1) / tail(r0, 1))
    for (L in LANDMARKS) {
        i <- max(which(tm$day_end <= L + 1e-9))
        est[paste0("rd_d", L)] <- unname(r1[i] - r0[i])
        est[paste0("rr_d", L)] <- unname(r1[i] / r0[i])
    }
    ## [4] stima descrittiva a fine grace period. NON e' un controllo di
    ## falsificazione: il clone "no early" viene censurato al momento
    ## dell'inizio dell'infusione, quindi le due strategie divergono gia'
    ## dentro la finestra e il contrasto non e' nullo per costruzione.
    ig <- max(which(tm$day_end <= GRACE_H/24 + 1e-9))
    est["rd_end_of_grace"] <- unname(r1[ig] - r0[ig])
    list(est = est, curve = data.frame(day = tm$day_end, risk0 = r0, risk1 = r1, rd = r1 - r0))
}

## ---- [3] bilanciamento sotto i pesi finali -----------------------------------
wmean <- function(x, w) { ok <- !is.na(x); sum(w[ok]*x[ok]) / sum(w[ok]) }
wsd   <- function(x, w) { ok <- !is.na(x); x <- x[ok]; w <- w[ok]
    m <- sum(w*x)/sum(w); sqrt(sum(w*(x-m)^2) / (sum(w) - sum(w^2)/sum(w))) }
smd_c <- function(x, g, w) (wmean(x[g==1], w[g==1]) - wmean(x[g==0], w[g==0])) /
    sqrt((wsd(x[g==1], w[g==1])^2 + wsd(x[g==0], w[g==0])^2) / 2)
smd_b <- function(x, g, w) { p1 <- wmean(x[g==1], w[g==1]); p0 <- wmean(x[g==0], w[g==0])
    (p1 - p0) / sqrt((p1*(1-p1) + p0*(1-p0)) / 2) }

## un clone-persona per braccio nel risk set post-grace, con il peso finale
final_balance <- function(pt, dat, cont_vars, bin_vars) {
    one <- pt[pt$in_grace == 0, ]; one <- one[!duplicated(one$key), ]
    g <- one$arm; w <- one$w; i <- one$row
    rbind(
        do.call(rbind, lapply(cont_vars, function(v)
            data.frame(variable = v, type = "cont",
                       smd_final = round(smd_c(dat[[v]][i], g, w), 3)))),
        do.call(rbind, lapply(bin_vars, function(v)
            data.frame(variable = v, type = "bin",
                       smd_final = round(smd_b(dat[[v]][i], g, w), 3)))))
}

weight_diag <- function(pt, g) {
    sc <- sat_curve(pt, g)
    ti <- attr(pt, "trunc")
    do.call(rbind, lapply(0:1, function(a) {
        s   <- pt$arm == a
        one <- pt[s & pt$in_grace == 0, ]; one <- one[!duplicated(one$key), ]
        data.frame(arm = a,
                   persons_after_grace = nrow(one),
                   censored_clones     = sum(pt$cens[s]),
                   events_unweighted   = sum(pt$event[s]),
                   ess_after_grace     = round(ess(one$w), 1),
                   weighted_risk_28d   = round(sc[[paste0("risk", a)]][nrow(sc)], 4),
                   w_q50 = round(median(one$w), 3),
                   w_q99 = round(unname(quantile(one$w, .99)), 3),
                   w_max = round(max(one$w), 3),
                   w_cens_cap = if (is.null(ti)) NA_real_ else ti$w_cens_cap[ti$arm == a],
                   pct_rows_truncated = if (is.null(ti)) NA_real_ else ti$pct_rows_truncated[ti$arm == a])
    }))
}

run_once <- function(dat, grace_h = GRACE_H, base_w = "overlap", diag = FALSE,
                     w_trunc = W_TRUNC_Q) {
    g  <- make_grid(grace_h)
    dat$treat <- as.integer(dat$h_init <= grace_h)   # etichetta per il PS soltanto
    pt <- build_weights(dat, build_pt(dat, grace_h, g), grace_h, base_w, w_trunc)
    out <- estimate(pt, g)
    if (diag) {
        out$diag <- weight_diag(pt, g)
        out$sat  <- sat_curve(pt, g)
        out$pt   <- pt
        out$n_at_risk_ps <- attr(pt, "n_at_risk_ps")
    }
    out
}

## =============================================================================
## 5. Imputazione multipla
## =============================================================================
make_imputations <- function(d, seed) {
    d$treat <- as.integer(d$h_init <= GRACE_H)
    iv <- unique(c(PS_VARS, "treat", "death_28")); iv <- iv[iv %in% names(d)]
    me <- make.method(d[, iv])
    me[intersect(c("treat", "death_28"), iv)] <- ""
    for (v in iv) if (!anyNA(d[[v]])) me[v] <- ""
    if (all(me == "")) return(list(data = replicate(M_IMP, d, simplify = FALSE), method = me))
    mi <- mice(d[, iv], m = M_IMP, method = me,
               predictorMatrix = make.predictorMatrix(d[, iv]),
               printFlag = FALSE, seed = seed)
    list(data = lapply(seq_len(M_IMP), function(m) {
        bm <- d; cm <- complete(mi, m)
        for (v in names(me)[me != ""]) bm[[v]] <- cm[[v]]
        bm
    }), method = me)
}

cat("\nProporzione mancante dopo il filtro di range:\n")
print(round(sapply(PS_VARS, function(v) mean(is.na(b0[[v]]))), 3))

MI      <- make_imputations(b0, SEED)
imputed <- MI$data
meth    <- MI$method
cat("Imputate:", paste(names(meth)[meth != ""], collapse = ", "), "\n")
b0$treat <- as.integer(b0$h_init <= GRACE_H)

## =============================================================================
## 6. Stima primaria: MI + bootstrap (ri-stima di PS e pesi), Rubin
## =============================================================================
t_start <- Sys.time()
res <- lapply(seq_len(M_IMP), function(m) {
    bm  <- imputed[[m]]; n <- nrow(bm)
    fit <- run_once(bm, diag = TRUE)
    bt  <- vector("list", B_BOOT)
    for (b in seq_len(B_BOOT)) {
        s <- sample.int(n, n, replace = TRUE)          # ricampionamento per paziente
        bt[[b]] <- tryCatch(run_once(bm[s, ]), error = function(e) NULL)
    }
    ok <- !vapply(bt, is.null, logical(1))
    bt_est <- do.call(rbind, lapply(bt[ok], `[[`, "est"))
    curve_var <- as.data.frame(lapply(c("risk0", "risk1", "rd"), function(k)
        apply(sapply(bt[ok], function(x) x$curve[[k]]), 1, var)))
    names(curve_var) <- c("risk0", "risk1", "rd")
    cat("imputazione", m, "/", M_IMP, "  RR =", round(fit$est[["rr"]], 4),
        "  bootstrap validi:", sum(ok), "\n")
    list(est = fit$est, var = apply(bt_est, 2, var), n_ok = sum(ok),
         curve = fit$curve, curve_var = curve_var, sat = fit$sat,
         n_at_risk_ps = fit$n_at_risk_ps,
         pt = if (m == 1) fit$pt else NULL,
         diag = cbind(imputation = m, fit$diag))
})
cat("Tempo:", round(difftime(Sys.time(), t_start, units = "mins"), 1), "minuti\n")

pool <- function(theta, within) {
    M <- length(theta); qbar <- mean(theta); ubar <- mean(within)
    bv <- if (M > 1) var(theta) else 0
    tot <- ubar + (1 + 1/M) * bv
    df  <- if (bv > 0) (M - 1) * (1 + ubar / ((1 + 1/M) * bv))^2 else Inf
    cr  <- qt(0.975, df = max(df, 1))
    c(est = qbar, lo = qbar - cr*sqrt(tot), hi = qbar + cr*sqrt(tot),
      fmi = ((1 + 1/M) * bv) / tot)
}
pool_log <- function(theta, within) {
    l <- pool(log(theta), within / theta^2)
    c(est = exp(l[["est"]]), lo = exp(l[["lo"]]), hi = exp(l[["hi"]]), fmi = l[["fmi"]])
}
get  <- function(nm) sapply(res, function(x) x$est[[nm]])
getv <- function(nm) sapply(res, function(x) x$var[[nm]])

rows <- list(rd_28 = pool(get("rd"), getv("rd")),
             rr_28 = pool_log(get("rr"), getv("rr")))
for (L in LANDMARKS) {
    rows[[paste0("rd_d", L)]] <- pool(get(paste0("rd_d", L)), getv(paste0("rd_d", L)))
    rows[[paste0("rr_d", L)]] <- pool_log(get(paste0("rr_d", L)), getv(paste0("rr_d", L)))
}
rows[["rd_end_of_grace"]] <- pool(get("rd_end_of_grace"), getv("rd_end_of_grace"))
primary <- do.call(rbind, lapply(names(rows), function(k)
    data.frame(estimate = k, value = rows[[k]][["est"]],
               lo = rows[[k]][["lo"]], hi = rows[[k]][["hi"]], fmi = rows[[k]][["fmi"]])))
risk0 <- mean(get("risk0")); risk1 <- mean(get("risk1"))

cat("\n=== PRIMARIA ===\n")
cat("Rischio a", HORIZON_D, "giorni: early", round(risk1, 4),
    "| no early", round(risk0, 4), "\n")
print(primary, row.names = FALSE)

gr <- primary[primary$estimate == "rd_end_of_grace", ]
cat("\n[Descrittiva] RD a fine grace (", GRACE_H, "h): ",
    sprintf("%+.4f (%+.4f a %+.4f)", gr$value, gr$lo, gr$hi),
    "\n  Non e' un contrasto nullo per costruzione: il clone no-early e'\n",
    "  censurato all'inizio dell'infusione, quindi le strategie divergono\n",
    "  gia' dentro la finestra.\n", sep = "")

ev <- function(rr) { r <- if (rr >= 1) rr else 1/rr; r + sqrt(r*(r - 1)) }
rr28 <- primary[primary$estimate == "rr_28", ]
evals <- c(ev(rr28$value), ev(rr28$lo))     # su valori NON arrotondati
cat("E-value stima:", round(evals[1], 2),
    "| E-value limite piu' vicino al nullo:", round(evals[2], 2), "\n")
write.csv(primary, file.path(out_dir, "primary_estimates.csv"), row.names = FALSE)

diag_all <- do.call(rbind, lapply(res, `[[`, "diag"))
write.csv(diag_all, file.path(out_dir, "weight_diagnostics.csv"), row.names = FALSE)
## [18] le colonne della troncatura sono tutte NA quando W_TRUNC_Q = 0, e
## aggregate() con la formula usa na.omit: senza questo filtro non resta
## nessuna riga da aggregare.
dcols <- setdiff(names(diag_all), "imputation")
dcols <- dcols[vapply(diag_all[dcols], function(x) !all(is.na(x)), TRUE)]
diag_mean <- aggregate(. ~ arm, data = diag_all[, dcols],
                       FUN = function(x) mean(x, na.rm = TRUE), na.action = na.pass)
cat("\n=== Diagnostiche dei pesi finali (media sulle imputazioni) ===\n")
print(round(diag_mean, 2), row.names = FALSE)

pool_curve <- function(k) {
    th <- sapply(res, function(x) x$curve[[k]]); wv <- sapply(res, function(x) x$curve_var[[k]])
    qbar <- rowMeans(th); tot <- rowMeans(wv) + (1 + 1/M_IMP) * apply(th, 1, var)
    data.frame(est = qbar, lo = qbar - 1.96*sqrt(tot), hi = qbar + 1.96*sqrt(tot))
}
cv <- data.frame(day = res[[1]]$curve$day)
for (k in c("risk0", "risk1", "rd")) {
    p <- pool_curve(k); cv[[k]] <- p$est; cv[[paste0(k, "_lo")]] <- p$lo; cv[[paste0(k, "_hi")]] <- p$hi
}
cv <- rbind(setNames(as.data.frame(matrix(0, 1, ncol(cv))), names(cv)), cv)
write.csv(cv, file.path(out_dir, "cumulative_incidence.csv"), row.names = FALSE)

## calibrazione: modello contro stimatore saturo
sat_bar <- Reduce(`+`, lapply(res, function(x) x$sat[, c("risk0", "risk1", "rd")])) / M_IMP
sat_bar$day <- res[[1]]$sat$day
cal <- data.frame(day = sat_bar$day,
                  risk0_model = cv$risk0[-1], risk0_saturated = sat_bar$risk0,
                  risk1_model = cv$risk1[-1], risk1_saturated = sat_bar$risk1)
cal$diff0 <- cal$risk0_model - cal$risk0_saturated
cal$diff1 <- cal$risk1_model - cal$risk1_saturated
write.csv(cal, file.path(chk_dir, "calibration_spline_vs_saturated.csv"), row.names = FALSE)
cat("\n=== Calibrazione spline vs saturo ===\n")
cat("Scarto assoluto massimo 0-28 g:", round(max(abs(c(cal$diff0, cal$diff1))), 4),
    "| oltre le prime 24 h:",
    round(max(abs(c(cal$diff0, cal$diff1)[rep(cal$day > 1, 2)])), 4), "\n")
cat("RR a 28 g — modello:", round(risk1 / risk0, 3),
    "| saturo:", round(tail(sat_bar$risk1, 1) / tail(sat_bar$risk0, 1), 3), "\n")
cat("NB: i due stimatori condividono gli stessi tempi di morte, quindi l'accordo\n",
    "non dice nulla sull'incertezza nella collocazione dei decessi.\n", sep = "")

## =============================================================================
## 7. Sensibilita'
##    [1] include la collocazione dei decessi noti come sola data (12 e 18 h;
##        il primario usa la fine del giorno, 24 h)
## =============================================================================
sens_row <- function(label, data_list = imputed, ...) {
    e <- lapply(data_list, function(d) run_once(d, ...)$est)
    data.frame(analysis = label,
               rr = exp(mean(log(sapply(e, `[[`, "rr")))),
               rd = mean(sapply(e, `[[`, "rd")))
}
imputed_raw <- make_imputations(b0_raw, SEED)$data
sens <- rbind(
    data.frame(analysis = "Primary (overlap, grace 12 h)",
               rr = rr28$value, rd = primary$value[primary$estimate == "rd_28"]),
    sens_row("No overlap tilt (full eligible cohort, different target population)",
             base_w = "none"),
    sens_row("No covariate range filtering", data_list = imputed_raw),
    sens_row(sprintf("Censoring weights truncated at the %.0fth percentile within arm",
                     100*W_TRUNC_SENS), w_trunc = W_TRUNC_SENS))
for (gh in GRACE_SENS)
    sens <- rbind(sens, sens_row(sprintf("Grace period %d h", gh), grace_h = gh))
if (n_date_only > 0) for (off in DOD_OFFSET_SENS) {
    dl <- lapply(imputed, set_death_time, offset_h = off)
    nc <- attr(dl[[1]]$h_to_death, "n_clipped")
    sens <- rbind(sens, sens_row(
        sprintf("Date-only deaths placed at %d h of the day of death (%d events anchored)",
                off, nc),
        data_list = dl))
}
sens$rr <- round(sens$rr, 3); sens$rd <- round(sens$rd, 4)
cat("\n=== SENSIBILITA' ===\n"); print(sens, row.names = FALSE)
write.csv(sens, file.path(out_dir, "sensitivity.csv"), row.names = FALSE)

## =============================================================================
## 8. Tabella 1, bilanciamento basale e bilanciamento finale
## =============================================================================
## Il punteggio usato per Tabella 1 e Figura 2 e' quello del tilt, stimato
## sull'intera coorte analitica: e' quello che definisce la pseudo-popolazione
## descritta nella tabella.
ps_bar <- rowMeans(sapply(imputed, function(d) {
    d$treat <- as.integer(d$h_init <= GRACE_H)
    f <- quiet_glm(glm(PS_FORM, family = binomial, data = d))
    pmin(pmax(as.numeric(fitted(f)), TRUNC_P), 1 - TRUNC_P)
}))
b0$ps <- ps_bar
b0$w  <- ifelse(b0$treat == 1, 1 - ps_bar, ps_bar)

LAB_C <- c(age = "Age, years", pf = "PaO2/FiO2 ratio, mm Hg", peep = "PEEP, cm H2O",
           fio2 = "FiO2, %", pao2 = "PaO2, mm Hg", lactate = "Lactate, mmol/L",
           creat = "Creatinine, mg/dL", plt = "Platelet count, K/uL",
           plateau = "Plateau pressure, cm H2O", driving_pressure = "Driving pressure, cm H2O",
           vt = "Tidal volume, mL", vt_pbw = "Tidal volume, mL/kg PBW",
           ne_equiv = "Norepinephrine equivalents, mcg/kg/min")
LAB_B <- c(female = "Female sex", vaso_active = "Vasopressor infusion at baseline",
           inotrope_active = "Inotrope infusion at baseline")
b0$female <- as.integer(b0$gender == "F")
for (lv in levels(b0$anchor_year_group)) {
    nm <- paste0("ayg_", gsub("[^0-9]", "", lv))
    b0[[nm]] <- as.integer(b0$anchor_year_group == lv)
    LAB_B[nm] <- paste("Anchor year group", lv)
}
LAB_C <- LAB_C[names(LAB_C) %in% names(b0)]; LAB_B <- LAB_B[names(LAB_B) %in% names(b0)]

tab1_row <- function(v, bin) {
    x <- b0[[v]]; g <- b0$treat; w1 <- rep(1, nrow(b0))
    f <- if (bin) smd_b else smd_c
    fmt <- function(w, grp) {
        m <- wmean(x[g == grp], w[g == grp])
        if (bin) sprintf("%.1f", 100*m) else sprintf("%.1f (%.1f)", m, wsd(x[g == grp], w[g == grp]))
    }
    data.frame(variable = unname(c(LAB_C, LAB_B)[v]),
               summary = if (bin) "%" else "mean (SD)",
               in_ps = v %in% PS_VARS |
                   (v == "female" && "gender" %in% PS_VARS) |
                   (grepl("^ayg_", v) && "anchor_year_group" %in% PS_VARS),
               n_missing = sum(is.na(x)),
               unwt_early = fmt(w1, 1), unwt_no_early = fmt(w1, 0),
               unwt_smd = round(abs(f(x, g, w1)), 3),
               wt_early = fmt(b0$w, 1), wt_no_early = fmt(b0$w, 0),
               wt_smd = round(abs(f(x, g, b0$w)), 3))
}
tab1 <- rbind(do.call(rbind, lapply(names(LAB_C), tab1_row, bin = FALSE)),
              do.call(rbind, lapply(names(LAB_B), tab1_row, bin = TRUE)))
tab1 <- rbind(tab1, data.frame(variable = "N (effective sample size, baseline overlap weights)",
                               summary = "n", in_ps = NA, n_missing = 0,
                               unwt_early = sum(b0$treat == 1), unwt_no_early = sum(b0$treat == 0),
                               unwt_smd = NA,
                               wt_early = round(ess(b0$w[b0$treat == 1]), 1),
                               wt_no_early = round(ess(b0$w[b0$treat == 0]), 1), wt_smd = NA))
cat("\n=== TABLE 1 ===\n"); print(tab1, row.names = FALSE)
write.csv(tab1, file.path(out_dir, "table1.csv"), row.names = FALSE)
ess_base <- c(early = ess(b0$w[b0$treat == 1]), no_early = ess(b0$w[b0$treat == 0]))
max_smd_ps <- max(tab1$wt_smd[which(tab1$in_ps)], na.rm = TRUE)
cat("|SMD| basale ponderata massima fra le covariate del PS:", round(max_smd_ps, 3), "\n")

## ---- [7] gravita' a T0 per finestra di inizio -------------------------------
## Se la gravita' basale e' simile fra gli strati ma la mortalita' cresce con
## il ritardo, il segnale non viene da cio' che era gia' vero a time zero.
win <- cut(b0$h_init, breaks = c(-Inf, 6, 12, 24, Inf),
           labels = c("0-6 h", "6-12 h", "12-24 h", ">24 h or never"), right = TRUE)
sev_vars <- intersect(c("pf", "peep", "plateau", "driving_pressure", "fio2",
                        "lactate", "ne_equiv", "age"), names(b0))
init_win <- data.frame(
    window = levels(win),
    n = as.integer(table(win)),
    pct = round(100 * as.numeric(table(win)) / nrow(b0), 1),
    death_28d_pct = round(100 * tapply(
        as.integer(!is.na(b0$h_to_death) & b0$h_to_death <= HORIZON_D*24),
        win, mean), 1))
for (v in sev_vars)
    init_win[[v]] <- round(tapply(b0[[v]], win, mean, na.rm = TRUE), 1)
cat("\n=== Gravita' a time zero e mortalita' per finestra di inizio ===\n")
print(init_win, row.names = FALSE)
write.csv(init_win, file.path(chk_dir, "severity_by_initiation_window.csv"),
          row.names = FALSE)

## [3] bilanciamento sotto i pesi finali, risk set post-grace, imputazione 1
b1 <- imputed[[1]]; b1$female <- as.integer(b1$gender == "F")
for (lv in levels(b0$anchor_year_group))
    b1[[paste0("ayg_", gsub("[^0-9]", "", lv))]] <- as.integer(b1$anchor_year_group == lv)
bal_fin <- final_balance(res[[1]]$pt, b1, names(LAB_C), names(LAB_B))
bal_fin$variable  <- unname(c(LAB_C, LAB_B)[bal_fin$variable])
bal_fin$smd_final <- abs(bal_fin$smd_final)
bal_fin <- merge(bal_fin, data.frame(variable = tab1$variable, smd_baseline = tab1$wt_smd),
                 by = "variable", all.x = TRUE)
cat("\n=== Bilanciamento sotto i pesi FINALI (risk set post-grace) ===\n")
print(bal_fin[order(-bal_fin$smd_final), c("variable", "smd_baseline", "smd_final")],
      row.names = FALSE)
write.csv(bal_fin, file.path(out_dir, "balance_final_weights.csv"), row.names = FALSE)
max_smd_fin <- max(bal_fin$smd_final, na.rm = TRUE)

nl <- cbind(b1$pf^2, b1$peep^2, b1$lactate^2, b1$pf*b1$vaso_active, b1$pf*b1$peep)
bal_nl <- data.frame(
    term = c("PaO2/FiO2 squared", "PEEP squared", "Lactate squared",
             "PaO2/FiO2 x vasopressor", "PaO2/FiO2 x PEEP"),
    smd_unweighted = round(abs(apply(nl, 2, smd_c, g = b0$treat, w = rep(1, nrow(b0)))), 3),
    smd_weighted   = round(abs(apply(nl, 2, smd_c, g = b0$treat, w = b0$w)), 3))
cat("\n=== Balance on terms not in the PS (baseline weights) ===\n")
print(bal_nl, row.names = FALSE)
write.csv(bal_nl, file.path(out_dir, "balance_nonlinear.csv"), row.names = FALSE)

## =============================================================================
## 9. Figure
## =============================================================================
ARM <- c("No early initiation", "Early initiation (within 12 h)")
COL <- setNames(c("#2166AC", "#B2182B"), ARM)
thm <- theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(), legend.position = "top")

stepify <- function(d) {
    n <- nrow(d)
    out <- d[c(rep(seq_len(n - 1), each = 2), n), , drop = FALSE]
    out$day <- d$day[c(1, rep(2:n, each = 2))]; out
}

keep <- tab1$summary != "n"
lp <- rbind(data.frame(variable = tab1$variable[keep], smd = tab1$unwt_smd[keep], phase = "Unweighted"),
            data.frame(variable = tab1$variable[keep], smd = tab1$wt_smd[keep],  phase = "Baseline overlap weighted"))
lp$variable <- factor(lp$variable, levels = tab1$variable[keep][order(tab1$unwt_smd[keep])])
lp$phase    <- factor(lp$phase, levels = c("Unweighted", "Baseline overlap weighted"))
ggsave(file.path(fig_dir, "fig1_love_baseline.png"),
       ggplot(lp, aes(smd, variable, colour = phase, shape = phase)) +
           geom_vline(xintercept = .1, linetype = "dashed", colour = "grey50") +
           geom_point(size = 2.4) +
           scale_colour_manual(values = c("grey35", "#B2182B")) +
           labs(x = "Absolute standardized mean difference", y = NULL,
                colour = NULL, shape = NULL) + thm,
       width = 6.5, height = 6.5, dpi = 300)

## [3] secondo love plot: pesi finali contro pesi basali
lp2 <- rbind(data.frame(variable = bal_fin$variable, smd = bal_fin$smd_baseline,
                        phase = "Baseline overlap weights"),
             data.frame(variable = bal_fin$variable, smd = bal_fin$smd_final,
                        phase = "Final clone-censor weights"))
lp2$variable <- factor(lp2$variable, levels = bal_fin$variable[order(bal_fin$smd_final)])
lp2$phase    <- factor(lp2$phase, levels = c("Baseline overlap weights", "Final clone-censor weights"))
ggsave(file.path(fig_dir, "figS_love_final.png"),
       ggplot(lp2, aes(smd, variable, colour = phase, shape = phase)) +
           geom_vline(xintercept = .1, linetype = "dashed", colour = "grey50") +
           geom_point(size = 2.4) +
           scale_colour_manual(values = c("grey35", "#2166AC")) +
           labs(x = "Absolute standardized mean difference", y = NULL,
                colour = NULL, shape = NULL,
                caption = "Final weights: one clone per person in the post-grace risk set.") + thm,
       width = 6.5, height = 6.5, dpi = 300)

dps <- data.frame(ps = b0$ps, arm = factor(b0$treat, 0:1, ARM))
ggsave(file.path(fig_dir, "fig2_ps.png"),
       ggplot(dps, aes(ps)) +
           geom_histogram(data = subset(dps, arm == ARM[2]),
                          aes(y = after_stat(density), fill = ARM[2]), bins = 40, alpha = .85) +
           geom_histogram(data = subset(dps, arm == ARM[1]),
                          aes(y = -after_stat(density), fill = ARM[1]), bins = 40, alpha = .85) +
           geom_hline(yintercept = 0, colour = "grey30") +
           scale_fill_manual(values = COL, breaks = ARM) +
           scale_y_continuous(labels = abs) +
           labs(x = "Propensity score", y = "Density", fill = NULL) + thm,
       width = 6.5, height = 4.5, dpi = 300)

curve <- rbind(
    stepify(data.frame(day = cv$day, risk = cv$risk0, lo = cv$risk0_lo, hi = cv$risk0_hi, arm = ARM[1])),
    stepify(data.frame(day = cv$day, risk = cv$risk1, lo = cv$risk1_lo, hi = cv$risk1_hi, arm = ARM[2])))
ggsave(file.path(fig_dir, "fig3_cif.png"),
       ggplot(curve, aes(day, risk, colour = arm, fill = arm)) +
           geom_ribbon(aes(ymin = pmax(lo, 0), ymax = hi), alpha = .18, colour = NA) +
           geom_line(linewidth = .9) +
           scale_colour_manual(values = COL) + scale_fill_manual(values = COL) +
           scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, NA)) +
           scale_x_continuous(breaks = seq(0, HORIZON_D, 7)) +
           labs(x = "Days since time zero", y = "Cumulative incidence of death (95% CI)",
                colour = NULL, fill = NULL) + thm,
       width = 6.5, height = 4.5, dpi = 300)

rdc <- stepify(data.frame(day = cv$day, rd = 100*cv$rd, lo = 100*cv$rd_lo, hi = 100*cv$rd_hi))
ggsave(file.path(fig_dir, "fig4_rd.png"),
       ggplot(rdc, aes(day, rd)) +
           geom_hline(yintercept = 0, colour = "grey40", linetype = "dashed") +
           geom_ribbon(aes(ymin = lo, ymax = hi), alpha = .18, fill = "#B2182B") +
           geom_line(linewidth = .9, colour = "#B2182B") +
           scale_x_continuous(breaks = seq(0, HORIZON_D, 7)) +
           labs(x = "Days since time zero",
                y = "Risk difference, percentage points (95% CI)") + thm,
       width = 6.5, height = 4, dpi = 300)

calp <- rbind(data.frame(day = sat_bar$day, sat = sat_bar$risk0, mod = cal$risk0_model, arm = ARM[1]),
              data.frame(day = sat_bar$day, sat = sat_bar$risk1, mod = cal$risk1_model, arm = ARM[2]))
ggsave(file.path(fig_dir, "figS_calibration.png"),
       ggplot(calp, aes(day, sat, colour = arm)) +
           geom_step(linewidth = .8) +
           geom_line(aes(y = mod), linetype = "22", linewidth = .8) +
           scale_colour_manual(values = COL) +
           scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
           scale_x_continuous(breaks = seq(0, HORIZON_D, 7)) +
           labs(x = "Days since time zero", y = "Cumulative incidence of death", colour = NULL,
                caption = paste("Solid: saturated weighted hazards.",   # [5]
                                "Dashed: pooled binomial model with complementary log-log link.")) +
           thm,
       width = 6.5, height = 4.5, dpi = 300)

## =============================================================================
## 10. Registro degli avvisi e manifest
## =============================================================================
wl <- WARN_LOG$msgs
writeLines(if (length(wl)) c(sprintf("%d avvisi non benigni:", length(wl)), names(table(wl)))
           else "Nessun avviso non benigno da glm (convergenza regolare in tutti i fit).",
           file.path(out_dir, "warnings_log.txt"))
cat("\nAvvisi non benigni registrati:", length(wl), "\n")

writeLines(c(
    paste("date:", Sys.time()),
    "script: 06_primary_ccw_final_v10.R",
    paste("cohort_file:", basename(COH)),
    paste("unit_of_analysis:", if (FIRST_STAY_ONLY) "first qualifying ICU stay per patient"
          else "ICU stay (patients may contribute more than one)"),
    paste("qualifying_stays:", n_stay_all, "| distinct patients:", n_subj_all,
          "| retained after restriction:", nrow(raw)),
    if (!is.null(TV_DIAG)) paste("time_varying_diagnostics:",
          paste(sprintf("%s: NA %.1f%%, differisce dal basale %.1f%%", TV_DIAG$variable,
                        100*TV_DIAG$pct_hours_na,
                        100*TV_DIAG$pct_hours_differs_from_baseline), collapse = "; "))
    else "time_varying_diagnostics: n/a",
    paste("censoring_model_covariates:",
          if (length(TV_USE)) paste("baseline plus time-varying:",
                                    paste(TV_USE, collapse = ", "))
          else "baseline only (hourly_grace.parquet not available)"),
    paste("seed:", SEED), paste("M:", M_IMP), paste("B:", B_BOOT),
    paste("bootstrap validi per imputazione:", paste(sapply(res, `[[`, "n_ok"), collapse = ", ")),
    paste("grace_h:", GRACE_H), paste("grace_sens_h:", paste(GRACE_SENS, collapse = ", ")),
    paste("horizon_d:", HORIZON_D), paste("landmarks_d:", paste(LANDMARKS, collapse = ", ")),
    paste("ps_trunc:", TRUNC_P),
    paste("censoring_weight_truncation:",
          if (W_TRUNC_Q > 0 && W_TRUNC_Q < 1) sprintf("quantile %.3f within arm", W_TRUNC_Q)
          else "none"),
    paste("ps_vars:", paste(PS_VARS, collapse = ", ")),
    paste("ps_tilt: fitted on the full analytic cohort (n =", nrow(b0), ")"),
    paste("ps_adherence: fitted in the risk set at the grace deadline (n =",
          res[[1]]$n_at_risk_ps, ")"),
    paste("imputed_vars:", paste(names(meth)[meth != ""], collapse = ", ")),
    paste("range_filtered_to_na:", paste(sprintf("%s=%d", impl$variable, impl$n_set_to_na),
                                         collapse = ", ")),
    paste("outcome: all-cause death;",
          if (HAS_DEATHTIME) "deathtime when available," else "deathtime NOT available,",
          "date-only deaths placed at", DOD_OFFSET_H, "h of the day of death;",
          "no censoring at discharge"),
    paste("deaths_date_only_all:", n_date_only,
          "| deaths_date_only_within_28d:", n_date_only_28,
          "| events anchored to first interval:", n_clipped),
    paste("deathtime_before_time_zero_set_to_missing:", n_dt_bad),
    paste("n_eligible:", flow$n[1]), paste("n_analytic:", nrow(b0)),
    paste("n_early:", sum(b0$treat)), paste("n_no_early:", sum(b0$treat == 0)),
    paste("ess_baseline_overlap: early", round(ess_base[["early"]], 1),
          "| no early", round(ess_base[["no_early"]], 1)),
    paste("ess_after_grace_final_weights: no early",
          round(diag_mean$ess_after_grace[diag_mean$arm == 0], 1),
          "| early", round(diag_mean$ess_after_grace[diag_mean$arm == 1], 1)),
    paste("max_abs_smd_baseline_ps_covariates:", round(max_smd_ps, 3)),
    paste("max_abs_smd_final_weights:", round(max_smd_fin, 3)),
    paste("risk_28d: early", round(risk1, 4), "| no early", round(risk0, 4)),
    paste("rr_28d:", round(rr28$value, 4), sprintf("(%.4f-%.4f)", rr28$lo, rr28$hi)),
    paste("rd_end_of_grace (descriptive, not a falsification test):",
          sprintf("%+.4f (%+.4f to %+.4f)", gr$value, gr$lo, gr$hi)),
    paste("calibration_max_abs_diff_after_24h:",
          round(max(abs(c(cal$diff0, cal$diff1)[rep(cal$day > 1, 2)])), 4)),
    paste("e_value:", round(evals[1], 2), "| e_value_ci:", round(evals[2], 2),
          "(computed from unrounded estimates)"),
    paste("non_benign_glm_warnings:", length(wl)),
    "", capture.output(sessionInfo())),
    file.path(out_dir, "run_manifest.txt"))

cat("\nTutto in:", out_dir, "\n")
