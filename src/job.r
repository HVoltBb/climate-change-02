"
Author/maintainer: Can Zhou [eidotog@gmail.com]
Date: Oct 4 2021
Version: 0.2
"

job <- function(i) {
    dyn.load(TMB::dynlib("v5_3"))
    set.seed(i)
    nlls <- matrix(NA, nrow = nrep, ncol = 10)
    indx <- 1:dim(tmpx)[1]

    for (repx in 1:nrep) {
        folds_n <- 10
        folds <- caret::createFolds(1:dim(tmpx)[1], k = folds_n)

        nll <- rep(NA, folds_n)

        for (fds in 1:folds_n) {
            test <- folds[fds]

            datx <- tmpx[indx[-unlist(test)], ]

            w <- table(datx$SpecimenID)
            dupes <- names(which(w > 1))
            singl <- names(which(w == 1))

            duplicity <- rep(0, dim(datx)[1])
            duplicity[which(datx$SpecimenID %in% singl)] <- 1:length(singl)

            for (i in 1:length(dupes)) {
                tmp <- which(datx$SpecimenID == dupes[i])
                duplicity[tmp] <- i + length(singl)
            }

            obj <- TMB::MakeADFun(
                data = list(
                    l1 = as.numeric(datx$ReLenCM),
                    l2 = as.numeric(datx$RcLenCM),
                    t1 = datx$REL_Date,
                    t2 = datx$REC_Date,
                    d1 = datx$d1,
                    d2 = datx$d2,
                    d1_d = datx$d1.d,
                    d2_d = datx$d2.d,
                    DUPE = duplicity - 1,
                    X_nao = X,
                    X_ao = X_a,
                    X_pna = X_p,
                    S = Matrix::.bdiag(list(S, S_a, S_p)),
                    Sdim = c(Sdim, Sdim_a, Sdim_p),
                    Pdim = dim(P_design)[1],
                    prediction_design_matrix_nao = P_design,
                    prediction_design_matrix_ao = P_design_a,
                    prediction_design_matrix_pna = P_design_p,
                    MNS = as.numeric(datx$RelMon) - 1,
                    MNS_map = nao.d$V2 - 1,
                    days = n_days.c / 30,
                    days_p = n_days,
                    reg = datx$reg,
                    p_dm_nao = p_dm_nao,
                    p_dm_ao = p_dm_ao,
                    p_dm_pna = p_dm_pna,
                    REFLEN = 80,
                    REFDAYS = 30,
                    REFMON = 7,
                    REFREG = 0,
                    STRTMON = 0,
                    PAR1 = PAR1, # 0: von Bertalanffy, 1: G. West, 2: Gompertz, 3: Logistic, 4: G-logistic, 4+: General von Bertalanffy
                    PAR2 = 1, # Deprecated. No effect.
                    PAR3 = .0, # Only considered when PAR1 > 4
                    PAR4 = 0, # break point between daily resolution and monthly resolution
                    PAR5 = 1, # 1: linf, 0: k [location of intrinsic effects]
                    PAR6 = 1 # 1: observation error,
                ),
                parameters = list(
                    c = c(298, 59, 47)[pickone],
                    logReg = 1,
                    reg_ = rep(0, 2),
                    k = c(0.008, 0.019, 0.025)[pickone],
                    logSig1 = 1.8,
                    gc = rep(0, Sdim),
                    logNaoc = -6.15,
                    gk = rep(0, Sdim),
                    logNaok = 0,
                    gc_ao = rep(0, Sdim_a),
                    logAoc = 1.1,
                    gk_ao = rep(0, Sdim_a),
                    logAok = 0,
                    gc_pna = rep(0, Sdim_p),
                    logPnac = 0.3,
                    gk_pna = rep(0, Sdim_p),
                    logPnak = 0,
                    e_o = rep(0, dim(datx)[1]),
                    logvInd = 3,
                    indl = rep(0, length(dupes) + length(singl)),
                    month_k = 0,
                    month_mg = 0,
                    lognu = 0
                ),
                map = list(
                    # Instructions:
                    # Comment OFF the each line to turn the corresponding effect ON,
                    # either in combination or individually except the month effect
                    # and seasonal pattern, which at most one of them can be turned
                    # ON simultaneously.
                    ####################################################################
                    logNaoc = factor(NA), gc = factor(rep(NA, Sdim)), # non-linear effect on linf
                    logNaok = factor(NA), gk = factor(rep(NA, Sdim)), # non-linear effect on k

                    logAoc = factor(NA), gc_ao = factor(rep(NA, Sdim_a)), # non-linear effect on linf
                    logAok = factor(NA), gk_ao = factor(rep(NA, Sdim_a)), # non-linear effect on k

                    logPnac = factor(NA), gc_pna = factor(rep(NA, Sdim_p)), # non-linear effect on linf
                    logPnak = factor(NA), gk_pna = factor(rep(NA, Sdim_p)), # non-linear effect on k

                    #                e_o  = factor(rep(NA, dim(datx)[1])),

                    month_k = factor(NA), month_mg = factor(NA),
                    logReg = factor(NA), reg_ = factor(rep(NA, 2)),
                    indl = factor(rep(NA, length(dupes) + length(singl))), logvInd = factor(NA), # individual effect
                    lognu = factor(NA),
                    ####################################################################
                    c = factor(1)
                ), # Do not edit this line. Just a placeholder.
                random = c("gc", "gk", "e_o", "indl", "gc_ao", "gk_ao", "gc_pna", "gk_pna", "reg_"),
                DLL = "v5_3",
                silent = TRUE
            )

            fit <- nlminb(fit_$par, obj$fn, obj$gr)
            TMBhelper::TMBAIC(fit)
            rep <- TMB::sdreport(obj)

            datx <- tmpx[indx[unlist(test)], ]

            w <- table(datx$SpecimenID)
            dupes <- names(which(w > 1))
            singl <- names(which(w == 1))

            duplicity <- rep(0, dim(datx)[1])
            duplicity[which(datx$SpecimenID %in% singl)] <- 1:length(singl)

            for (i in 1:length(dupes)) {
                tmp <- which(datx$SpecimenID == dupes[i])
                duplicity[tmp] <- i + length(singl)
            }

            obj <- TMB::MakeADFun(
                data = list(
                    l1 = as.numeric(datx$ReLenCM),
                    l2 = as.numeric(datx$RcLenCM),
                    t1 = datx$REL_Date,
                    t2 = datx$REC_Date,
                    d1 = datx$d1,
                    d2 = datx$d2,
                    d1_d = datx$d1.d,
                    d2_d = datx$d2.d,
                    DUPE = duplicity - 1,
                    X_nao = X,
                    X_ao = X_a,
                    X_pna = X_p,
                    S = Matrix::.bdiag(list(S, S_a, S_p)),
                    Sdim = c(Sdim, Sdim_a, Sdim_p),
                    Pdim = dim(P_design)[1],
                    prediction_design_matrix_nao = P_design,
                    prediction_design_matrix_ao = P_design_a,
                    prediction_design_matrix_pna = P_design_p,
                    MNS = as.numeric(datx$RelMon) - 1,
                    MNS_map = nao.d$V2 - 1,
                    days = n_days.c / 30,
                    days_p = n_days,
                    reg = datx$reg,
                    p_dm_nao = p_dm_nao,
                    p_dm_ao = p_dm_ao,
                    p_dm_pna = p_dm_pna,
                    REFLEN = 80,
                    REFDAYS = 30,
                    REFMON = 7,
                    REFREG = 0,
                    STRTMON = 0,
                    PAR1 = PAR1, # 0: von Bertalanffy, 1: G. West, 2: Gompertz, 3: Logistic, 4: G-logistic, 4+: General von Bertalanffy
                    PAR2 = 1, # Deprecated. No effect.
                    PAR3 = .0, # Only considered when PAR1 > 4
                    PAR4 = 0, # break point between daily resolution and monthly resolution
                    PAR5 = 1, # 1: linf, 0: k [location of intrinsic effects]
                    PAR6 = 1 # 1: observation error,
                ),
                parameters = list(
                    c = c(256, 59, 47)[pickone],
                    logReg = 0.96,
                    reg_ = rep(0, 2),
                    k = c(0.01, 0.019, 0.025)[pickone],
                    logSig1 = 1.8,
                    gc = rep(0, Sdim),
                    logNaoc = -6.15,
                    gk = rep(0, Sdim),
                    logNaok = 0,
                    gc_ao = rep(0, Sdim_a),
                    logAoc = 1.1,
                    gk_ao = rep(0, Sdim_a),
                    logAok = 0,
                    gc_pna = rep(0, Sdim_p),
                    logPnac = 0.3,
                    gk_pna = rep(0, Sdim_p),
                    logPnak = 0,
                    e_o = rep(0, dim(datx)[1]),
                    logvInd = 2.3,
                    indl = rep(0, length(dupes) + length(singl)),
                    month_k = 0,
                    month_mg = 0,
                    lognu = 0
                ),
                map = list(
                    # Instructions:
                    # Comment OFF the each line to turn the corresponding effect ON,
                    # either in combination or individually except the month effect
                    # and seasonal pattern, which at most one of them can be turned
                    # ON simultaneously.
                    ####################################################################
                    logNaoc = factor(NA), gc = factor(rep(NA, Sdim)), # non-linear effect on linf
                    logNaok = factor(NA), gk = factor(rep(NA, Sdim)), # non-linear effect on k

                    logAoc = factor(NA), gc_ao = factor(rep(NA, Sdim_a)), # non-linear effect on linf
                    logAok = factor(NA), gk_ao = factor(rep(NA, Sdim_a)), # non-linear effect on k

                    logPnac = factor(NA), gc_pna = factor(rep(NA, Sdim_p)), # non-linear effect on linf
                    logPnak = factor(NA), gk_pna = factor(rep(NA, Sdim_p)), # non-linear effect on k

                    #                e_o  = factor(rep(NA, dim(datx)[1])),

                    month_k = factor(NA), month_mg = factor(NA),
                    logReg = factor(NA), reg_ = factor(rep(NA, 2)),
                    logvInd = factor(NA), indl = factor(rep(NA, length(dupes) + length(singl))), # individual effect
                    lognu = factor(NA),
                    ####################################################################
                    c = factor(1)
                ), # Do not edit this line. Just a placeholder.
                random = c("e_o", "indl"),
                DLL = "v5_3",
                silent = TRUE
            )

            nll[fds] <- obj$fn(c(
                c = fit$par["c"],
                k = fit$par["k"],
                #            gc = rep$par.random[names(rep$par.random)=="gc"],
                #            gc_ao = rep$par.random[names(rep$par.random)=="gc_ao"],
                #            gc_pna = rep$par.random[names(rep$par.random)=="gc_pna"],
                logSig1 = fit$par["logSig1"]
                #            month_k = fit$par['month_k'], month_mg = fit$par['month_mg'],
                #            logvInd = fit$par['logvInd'],
                #            reg_ = rep$par.random[names(rep$par.random)=='reg_']
                #            lognu = fit$par['lognu']
            ))
        }

        nlls[repx, ] <- nll
    }
    return(nlls)
}