Frailty.h<-function(formulaMain,censor,DataMain,RandDist="Normal",mord=0,dord=1,Maxiter=200,convergence=1e-7,contrasts=NULL){
    require(Matrix)
    require(numDeriv)
    mc <- match.call()
    fr <- HGLMFrames(mc, formulaMain, contrasts)
    namesX <- names(fr$fixef)
    namesY <- names(fr$mf)[1]
    FL <- HGLMFactorList(formulaMain, fr, 0L, 0L)
    namesRE <- FL$namesRE
    y <- matrix(fr$Y, length(fr$Y), 1)
    x <- fr$X
    z <- FL$Design
    n<-nrow(x)
    p<-ncol(x)
    nrand <- length(z)
    q <- rep(0, nrand)
    for (i in 1:nrand) q[i] <- dim(z[[i]])[2]
    del <-matrix(0,n,1)
    del[,1] <- censor
    SS <- FL$Subject
    res1<-FrailtyMakeData(y,x,del,z)
    y<-res1[1][[1]]
    x<-res1[2][[1]]
    del<-res1[3][[1]]
    z<-res1[4][[1]]
    Mi<-res1[5][[1]]
    idx2<-res1[6][[1]]
    t2<-res1[7][[1]]
    di<-res1[8][[1]]
    beta_h<-matrix(0,p,1)
    qcum <- cumsum(c(0, q))
    v_h<-matrix(0,qcum[nrand+1],1)
    alpha_h <- rep(0, nrand)
    for (i in 1:nrand) alpha_h[i] <- 0.1
    Max_iter<-Maxiter
    err<-1
    for ( i in 1:Max_iter) {
        if (err>=0.000001) {
        if (RandDist=="Normal") res2<-PNFrailty.h(x,z,y, del,Mi,idx2,t2, di, beta_h,v_h, alpha_h,mord,dord)
        if (RandDist=="Gamma") res2<-PGFrailty.h(x,z,y, del,Mi,idx2,t2, di, beta_h,v_h, alpha_h,mord,dord)
        alpha_h<-res2[13][[1]]
        alpha_h1<-res2[14][[1]]
        temp4<-sum(abs(alpha_h-alpha_h1))
        err<-temp4
        beta_h<-res2[11][[1]]
        v_h<-res2[12][[1]]
        alpha_h<-alpha_h1
        se_beta<-res2[20][[1]]
        print_i<-i
        print_err<-err
        names(print_i) <- "iteration : "
        print(print_i)
        names(print_err) <- "convergence : "
        print(print_err)
        }
    }
###############################################################
############# print estimates ###########################
###############################################################
    if (RandDist=="Gamma") print("Results from the gamma frailty model")
    if (RandDist=="Normal") print("Results from the log-normal frailty model")
    print(formulaMain)
    if (mord==0 && dord==1) print("Method : HL(0,1)")   
    if (mord==0 && dord==2) print("Method : HL(0,2)")     
    if (mord==1 && dord==1) print("Method : HL(1,1)")   
    if (mord==1 && dord==2) print("Method : HL(1,2)")   
    print("Estimates from the mean model")
    z_beta<-beta_h/se_beta
    pval <- 2 * pnorm(abs(z_beta), lower.tail = FALSE)
    beta_coeff<-cbind(beta_h,se_beta,z_beta)
    colnames(beta_coeff) <- c("Estimate", "Std. Error", "t-value")
    rownames(beta_coeff) <- namesX
    print(beta_coeff,4)
###############################################################
############# se for lambda         ###########################
###############################################################
    if (RandDist=="Normal") res3<-PNFrailty_SE.h(res2,nrand,q,qcum,dord)
    if (RandDist=="Gamma") res3<-PGFrailty_SE.h(res2,nrand,q,qcum,dord)
    print("Estimates from the dispersion model")
    se_alpha_h<-res3[1][[1]]
    hlike<--2*res3[2][[1]]
    p1<--2*res3[3][[1]]
    p2<--2*res3[4][[1]]
    p3<--2*res3[5][[1]]
    z_lam<-alpha_h/se_alpha_h
    lam_coeff<-cbind(alpha_h,se_alpha_h)
    colnames(lam_coeff) <- c("Estimate", "Std. Error")
    rownames(lam_coeff) <- namesRE
    print(lam_coeff,4)
###############################################################
############# -2*Likelihoods         ##########################
###############################################################
    if (mord==0 && dord==1) like_value<-cbind(hlike,p1)
    if (mord==0 && dord==1) colnames(like_value) <- c("-2*hp","-2*p_b,v(hp)")
    if (mord==0 && dord==2) like_value<-cbind(hlike,p3)
    if (mord==0 && dord==2) colnames(like_value) <- c("-2*hp","-2*s_b,v(hp)")
    if (mord==1 && dord==1) like_value<-cbind(hlike,p2,p1)
    if (mord==1 && dord==1) colnames(like_value) <- c("-2*hp","-2*p_v(hp)","-2*p_b,v(hp)")
    if (mord==1 && dord==2) like_value<-cbind(hlike,p2,p3)
    if (mord==1 && dord==2) colnames(like_value) <- c("-2*hp","-2*p_v(hp)","-2*s_b,v(hp)")
    print(like_value,5)
    res4<-list(res2,res3)
    return(res4)   
}

FrailtyMakeData<-function(y,x,del,z) {
    n<-nrow(x)
    p<-ncol(x)
    nrand <- length(z)
    q <- rep(0, nrand)
    for (i in 1:nrand) q[i] <- dim(z[[i]])[2]
    zz<-z[[1]]
    if (nrand>1) {
        index1<-nrand
        for (i in 2:index1) zz<-cbind(zz, z[[i]])
    }
    qcum <- cumsum(c(0, q))
    zzz<-matrix(0,n,qcum[nrand+1])
    zzz[1:n,1:qcum[nrand+1]]<-zz[1:n,1:qcum[nrand+1]]
    sort_data<-cbind(y,x,del,zzz)
    sort.res<-sort_data[order(sort_data[,1],na.last=NA),]
    y[1:n,1]<-sort.res[1:n,1]
    index1<-p+1
    x[1:n,1:p]<-sort.res[1:n,2:index1]
    index1<-index1+1
    del[1:n,1]<-sort.res[1:n,index1]
    for (i in 1:nrand) {
        index3<-p+2+qcum[i]+1
        index4<-p+2+qcum[i+1]
        z[[i]][1:n,1:q[i]]<-sort.res[1:n,index3:index4]
    }

t<-matrix(0,n,1)
xx<-matrix(0,n,p)
di<-matrix(0,n,1)
idx1<-0

for (i in 1:n) {
if (del[i,1]==1) {
   idx1<-idx1+1
   t[idx1,1]<-y[i,1]
   for (j in 1:p) {
    xx[idx1,j]<-x[i,j]
   }
}
}

t1<-t
for (i in 1:idx1) {
for (j in 1:idx1) {
if (t1[i,1]==t1[j,1]) {
   if(i != j) {
    t1[j,1]<-0
   }
}
}
}

t2<-matrix(0,idx1,1)
idx2<-0
for (i in 1:idx1) {
if (t1[i,1] != 0) {
 idx2<-idx2+1
 t2[idx2,1]<-t1[i,1]
}
}

di<-matrix(0,idx2,1)
si<-matrix(0,idx2,p)
for (i in 1:idx2) {
   di[i,1]<-0
   for (j in 1:idx1) {
    if(t2[i,1]==t[j,1]) {
     di[i,1]<-di[i,1]+1
     for(k in 1:p) {
      si[i,k]<-si[i,k]+xx[j,k]
     }
    }
   }
}

# Triangle Mat form
Mi<-matrix(0,n,idx2)
for (i in 1:n) {
  t0<-y[i,1]
 for (j in 1:idx2) {
 if (t2[j,1] <= t0) {
  Mi[i,j]=1
 } else Mi[i,j]=0
 }
}
  res<-list(y,x,del,z,Mi,idx2,t2, di)
  return(res)
}




FrailtyMakeData1<-function(y,x,del,z) {
n<-nrow(x)
p<-ncol(x)
nrand <- length(z)
q <- rep(0, nrand)
for (i in 1:nrand) q[i] <- dim(z[[i]])[2]

# handle distinct death time

index1<-n-1
temp2<-matrix(0,n,p)
for (i in 1:index1){
index2<-n-i
for (j in 1:index2) {
   if(y[j,1]>y[j+1,1]) {
    temp<-y[j,1] 
    y[j,1]<-y[j+1,1] 
    y[j+1,1]<-temp 

    temp1<-del[j,1] 
    del[j,1]<-del[j+1,1] 
    del[j+1,1]<-temp1 

    for (k in 1:p) temp2[j,k]<-x[j,k]
    for ( k in 1:p) x[j,k]<-x[j+1,k]
    for ( k in 1:p) x[j+1,k]<-temp2[j,k]

    for (l in 1:nrand) {
       temp3<-matrix(0,n,q[l])
       for (k in 1:q[l]) temp3[j,k]<-z[[l]][j,k]
       for (k in 1:q[l]) z[[l]][j,k]<-z[[l]][j+1,k]
       for (k in 1:q[l]) z[[l]][j+1,k]<-temp3[j,k]
   }
  }
}
}

t<-matrix(0,n,1)
xx<-matrix(0,n,p)
di<-matrix(0,n,1)
idx1<-0

for (i in 1:n) {
if (del[i,1]==1) {
   idx1<-idx1+1
   t[idx1,1]<-y[i,1]
   for (j in 1:p) {
    xx[idx1,j]<-x[i,j]
   }
}
}

t1<-t
for (i in 1:idx1) {
for (j in 1:idx1) {
if (t1[i,1]==t1[j,1]) {
   if(i != j) {
    t1[j,1]<-0
   }
}
}
}

t2<-matrix(0,idx1,1)
idx2<-0
for (i in 1:idx1) {
if (t1[i,1] != 0) {
 idx2<-idx2+1
 t2[idx2,1]<-t1[i,1]
}
}

di<-matrix(0,idx2,1)
si<-matrix(0,idx2,p)
for (i in 1:idx2) {
   di[i,1]<-0
   for (j in 1:idx1) {
    if(t2[i,1]==t[j,1]) {
     di[i,1]<-di[i,1]+1
     for(k in 1:p) {
      si[i,k]<-si[i,k]+xx[j,k]
     }
    }
   }
}

# Triangle Mat form
Mi<-matrix(0,n,idx2)
for (i in 1:n) {
  t0<-y[i,1]
 for (j in 1:idx2) {
 if (t2[j,1] <= t0) {
  Mi[i,j]=1
 } else Mi[i,j]=0
 }
}
  res<-list(y,x,del,z,Mi,idx2,t2, di)
  return(res)
}


HGLMFrames<-function (mc, formula, contrasts, vnms = character(0)) 
{
    mf <- mc
    m <- match(c("DataMain", "weights", "na.action", "offset"), 
        names(mf), 0)
    mf <- mf[c(1, m)]
    frame.form <- subbars(formula)
    if (length(vnms) > 0) 
        frame.form[[3]] <- substitute(foo + bar, list(foo = parse(text = paste(vnms, 
            collapse = " + "))[[1]], bar = frame.form[[3]]))
    fixed.form <- nobars(formula)
    if (inherits(fixed.form, "name")) 
        fixed.form <- substitute(foo ~ 1, list(foo = fixed.form))
    environment(fixed.form) <- environment(frame.form) <- environment(formula)
    mf$formula <- frame.form
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    names(mf)[2] <- "data"
    fe <- mf
    mf <- eval(mf, parent.frame(2))
    fe$formula <- fixed.form
    fe <- eval(fe, parent.frame(2))
    fe
    Y <- model.response(mf, "any")
    if (length(dim(Y)) == 1) {
        nm <- rownames(Y)
        dim(Y) <- NULL
        if (!is.null(nm)) 
            names(Y) <- nm
    }
    mt <- attr(fe, "terms")
    X <- if (!is.empty.model(mt)) 
        model.matrix(mt, mf, contrasts)
    else matrix(, NROW(Y), 0)
    storage.mode(X) <- "double"
    fixef <- numeric(ncol(X))
    names(fixef) <- colnames(X)
    dimnames(X) <- NULL
    wts <- model.weights(mf)
    if (is.null(wts)) 
        wts <- numeric(0)
    off <- model.offset(mf)
    if (is.null(off)) 
        off <- numeric(0)
    if (any(wts <= 0)) 
        stop(gettextf("negative weights or weights of zero are not allowed"))
    if (length(off) && length(off) != NROW(Y)) 
        stop(gettextf("number of offsets is %d should equal %d (number of observations)", 
            length(off), NROW(Y)))
    attr(mf, "terms") <- mt
    list(Y = Y, X = X, wts = as.double(wts), off = as.double(off), 
        mf = mf, fixef = fixef)
}

HGLMFactorList <- function (formula, fr, rmInt, drop) 
{
    mf <- fr$mf
    bars <- expandSlash(findbars(formula[[3]]))
    for (i in 1:length(bars)) {
        checkcorr <- findplus(bars[[i]])
        if (checkcorr == 1) 
            stop("Correlated random effects are not currently allowed in the HGLM routines")
        if (checkcorr == -1) 
            stop("You do not need to specify '-1' for no intercept it is done be default")
    }
    if (!length(bars)) 
        stop("No random effects terms specified in formula")
    names(bars) <- unlist(lapply(bars, function(x) deparse(x[[3]])))
    fl <- lapply(bars, function(x) {
        ff <- eval(substitute(as.factor(fac)[, drop = TRUE], 
            list(fac = x[[3]])), mf)
        im <- as(ff, "sparseMatrix")
        if (!isTRUE(validObject(im, test = TRUE))) 
            stop("invalid conditioning factor in random effect: ", 
                format(x[[3]]))
        if (is.name(x[[2]])) {
            tempexp <- paste("~", as.character(x[[2]]), "-1")
            tempexp <- as.formula(tempexp)[[2]]
        }
        else tempexp <- x[[2]]
        mm <- model.matrix(eval(substitute(~expr, list(expr = tempexp))), 
            mf)
        if (rmInt) {
            if (is.na(icol <- match("(Intercept)", colnames(mm)))) 
                break
            if (ncol(mm) < 2) 
                stop("lhs of a random-effects term cannot be an intercept only")
            mm <- mm[, -icol, drop = FALSE]
        }
        ans <- list(f = ff, A = do.call(rBind, lapply(seq_len(ncol(mm)), 
            function(j) im)), Zt = do.call(rBind, lapply(seq_len(ncol(mm)), 
            function(j) {
                im@x <- mm[, j]
                im
            })), ST = matrix(0, ncol(mm), ncol(mm), dimnames = list(colnames(mm), 
            colnames(mm))))
        if (drop) {
            ans$A@x <- rep(0, length(ans$A@x))
            ans$Zt <- drop0(ans$Zt)
        }
        ans
    })
    Design <- list(0)
    Subject <- list(0)
    for (i in 1:length(fl)) {
        Subject[[i]] <- as.factor(fl[[i]]$f)
        tempmat <- fl[[i]]$Zt
        tempmat <- as.matrix(t(tempmat))
        Design[[i]] <- tempmat
    }
    list(Design = Design, Subject = Subject, namesRE = names(bars))
}

PNFrailty.h<-function(x,z,y,del,Mi,idx2,t2,di,beta_h0,v_h0,alpha_h0,mord,dord){
    n<-nrow(x)
    p<-ncol(x)
    nrand <- length(z)
    q <- rep(0, nrand)
    for (i in 1:nrand) q[i] <- dim(z[[i]])[2]
    qcum <- cumsum(c(0, q))
    beta_h<-beta_h0
    v_h<-v_h0
    alpha_h<-alpha_h0
    
    zz<-z[[1]]
    if (nrand>1) {
        index1<-nrand
        for (i in 2:index1) zz<-cbind(zz, z[[i]])
    }
    z<-matrix(0,n,qcum[nrand+1])
    z[1:n,1:qcum[nrand+1]]<-zz[1:n,1:qcum[nrand+1]]
    muh<-x%*%beta_h0 + z%*%v_h0

    expeta<-exp(muh)
    cla0<-di/(t(Mi)%*%expeta)
    Wi<-diag(expeta[,1])
    Ai<-diag(cla0[,1])
    done<-matrix(1,idx2,1)

########################
    oq<-matrix(1,qcum[nrand+1],1)
    oq1<-matrix(1,qcum[nrand+1],1)
    clam0<-Mi%*%Ai%*%done
    
    for (i in 1:nrand) {
       index1<-qcum[i]+1
       oq1[index1:qcum[i+1]]<-alpha_h[i]
    } 
    D<-diag(oq1[,1])
    iD<-solve(D)
    dft1<-t(x)%*%(del-Wi%*%clam0)
######################## pv(hp) for beta ########################
    if (mord==1) {
    U <- iD
    Bi<-diag(clam0[,1])
    temp4<-cla0^2/di
    As<-diag(temp4[,1])
    mat<-(Wi%*%Bi)-(Wi%*%Mi)%*%As%*%(t(Mi)%*%Wi)
    dinv<-solve(t(z)%*%mat%*%z+U)
    mu0<-exp(x%*%beta_h0 + z%*%v_h0)
    mu<-exp(x%*%beta_h0 + z%*%v_h0)*clam0
    dcla0be<--As%*%(t(Mi)%*%Wi)
    dcla0b<-matrix(0,idx2,p)
    dv_db<-matrix(0,qcum[nrand+1],p)
    xz<-matrix(0,n,p)
    dmu0<-matrix(0,n,p)
    dw_db1<-matrix(0,n,p)
    xk<-matrix(0,n,1)
    ad1<-matrix(0,p,1)
    for (k in 1:p) {
       xk[,1]<-x[,k]
       dv_db[,k] <--dinv%*%(t(z)%*%mat%*%xk)
       xz[,k]<-xk+z%*%(dv_db[,k])
       dcla0b[,k]<-dcla0be%*%(xz[,k])
       dc<-Mi%*%diag(dcla0b[,k])%*%done
       dmu0[,k]<-mu0*(xz[,k])
       dw_db1[,k]<- dc*mu0 + mu*xz[,k]
       temp4<-(2*cla0/di)*(dcla0b[,k])
       dw_db2<-((diag(dmu0[,k]))%*%Mi%*%As%*%t(Mi)%*%Wi)+(Wi%*%Mi%*%diag(temp4[,1])%*%t(Mi)%*%Wi)+(Wi%*%Mi%*%As%*%t(Mi)%*%(diag(dmu0[,k]))) 
       dw_db<-diag(dw_db1[,k])-dw_db2
       ad1[k,1]<-sum(diag(dinv%*%(t(z)%*%dw_db%*%z))) 
    }
    dft1<-dft1-0.5*ad1
   }
#########################################################
    dft2<-t(z)%*%(del-Wi%*%clam0)-(iD%*%v_h0)
    dft<-rbind(dft1,dft2)
    Bi<-diag(clam0[,1])
    temp4<-cla0^2/di
    Adi<-diag(temp4[,1])
    As<-Adi
    mat<-(Wi%*%Bi)-(Wi%*%Mi)%*%Adi%*%(t(Mi)%*%Wi)
    U<-iD
    H <- rbind(cbind(t(x)%*%mat%*%x, t(x)%*%mat%*%z), cbind(t(z)%*%mat%*%x, t(z)%*%mat%*%z+U))
    Hinv<-solve(H)
    be_h0<- rbind(beta_h0, v_h0)
    be_h <- be_h0 + (Hinv%*%dft)
    beta_h[1:p,1]<-be_h[1:p,1]
    se_beta_h<-matrix(0,p,1)
    for (i in 1:p) se_beta_h[i,1]<-sqrt(Hinv[i,i])
    index2<-qcum[nrand+1]
    index3<-p+1
    index4<-p+qcum[nrand+1]
    v_h[1:index2,1]<-be_h[index3:index4,1]
################################################
    for (i in 1:nrand) {
    if (dord==0) {
        index1<-p+qcum[i]+1
        index2<-p+qcum[i+1]
        gamma<-sum(diag(Hinv[index1:index2,index1:index2]))/alpha_h[i]
        index1<-qcum[i]+1
        index2<-qcum[i+1]
        alpha_h[i]<-sum(v_h[index1:index2,1]^2)/(q[i]-gamma)
    }
    if (dord==1 | dord==2) {
        H22<-solve(t(z)%*%mat%*%z+U)
        ial1<-1/alpha_h[i]
        iA<-iD
        C<-matrix(0,qcum[nrand+1],qcum[nrand+1])
        index1<-qcum[i]+1
        index2<-qcum[i+1]
        for (j in index1:index2) C[j,j]<-1
        iB1<-iA%*%C%*%iA
        c_vh1<-iB1%*%v_h
        dv1<-H22%*%c_vh1
        dexpeta1<-expeta*(z%*%dv1)
        dcla01<--(di/((t(Mi)%*%expeta)^2))*(t(Mi)%*%dexpeta1)
        dWi1<-diag(dexpeta1[,1])
        dAi1<-diag(dcla01[,1])
        temp4<-Mi%*%dAi1%*%done
        dBi1<-diag(temp4[,1])
        dvec1<-2*(cla0*dcla01)
        temp4<-dvec1/di
        dAs1<-diag(temp4[,1])
        dmat1<-(dWi1%*%Bi)+(Wi%*%dBi1)-(dWi1%*%Mi%*%As%*%t(Mi)%*%Wi)-(Wi%*%Mi%*%dAs1%*%t(Mi)%*%Wi)-(Wi%*%Mi%*%As%*%t(Mi)%*%dWi1)
        dia1<--iB1
        Hd1 <- rbind(cbind(t(x)%*%dmat1%*%x,t(x)%*%dmat1%*%z),cbind(t(z)%*%dmat1%*%x,t(z)%*%dmat1%*%z+ dia1))
        gamma1<--alpha_h[i]*sum(diag((Hinv%*%Hd1)))
        if (dord==2) {
        expeta <- exp(x%*%beta_h + z%*%v_h)
        ial1<-1/alpha_h[i]
        muu<-exp(x%*%beta_h)*clam0
        zmuu<-t(z)%*%muu
        u_h<-exp(v_h)
        ude1<-u_h*dv1
        aa1<-(zmuu*u_h)+ial1
        bb1<-(zmuu*ude1)-(ial1^2)
        term11<-((aa1*zmuu*ude1)-(2*zmuu*u_h*bb1))/(aa1^3)
        term21<-((2*aa1*zmuu*zmuu*u_h*ude1)-(3*((zmuu*u_h)^2)*bb1))/(aa1^4)
        term1<-(3*term11)-(5*term21)
        SS1<-diag(term1[,1])
        gamma21<--(alpha_h[i]/12)*sum(diag(SS1))
        }
        if (dord==1) {
            gamma21<-0
        }
        k21<- q[i]- gamma1- gamma21-sum(v_h[index1:index2,1]^2)/(alpha_h[i])
        alpha_h[i]<-sum(v_h[index1:index2,1]^2)/(q[i]-gamma1-gamma21)
    }
    }
    res<-list(x,z,y,del,Mi,idx2,t2,di,beta_h0,v_h0,beta_h,v_h,alpha_h0,alpha_h,dft,Hinv,clam0,H,mat,se_beta_h,U)
    return(res)
}

PGFrailty.h<-function(x,z,y,del,Mi,idx2,t2,di,beta_h0,v_h0,alpha_h0,mord,dord){
    n<-nrow(x)
    p<-ncol(x)
    nrand <- length(z)
    q <- rep(0, nrand)
    for (i in 1:nrand) q[i] <- dim(z[[i]])[2]
    qcum <- cumsum(c(0, q))
    beta_h<-beta_h0
    v_h<-v_h0
    alpha_h<-alpha_h0
    
    zz<-z[[1]]
    if (nrand>1) {
        index1<-nrand
        for (i in 2:index1) zz<-cbind(zz, z[[i]])
    }
    z<-matrix(0,n,qcum[nrand+1])
    z[1:n,1:qcum[nrand+1]]<-zz[1:n,1:qcum[nrand+1]]
    muh<-x%*%beta_h0 + z%*%v_h0
    expeta<-exp(muh)
    cla0<-di/(t(Mi)%*%expeta)
    Wi<-diag(expeta[,1])
    Ai<-diag(cla0[,1])
    done<-matrix(1,idx2,1)

########################
    oq<-matrix(1,qcum[nrand+1],1)
    oq1<-matrix(1,qcum[nrand+1],1)
    clam0<-Mi%*%Ai%*%done
    
    for (i in 1:nrand) {
       index1<-qcum[i]+1
       oq1[index1:qcum[i+1]]<-alpha_h[i]
    } 
    D<-diag(oq1[,1])
    iD<-solve(D)
    iu_h0<-exp(v_h0) ## gamma frailty
    U<-iD%*%diag(iu_h0[,1]) ## gamma frailty
    dft1<-t(x)%*%(del-Wi%*%clam0)
######################## pv(hp) for beta ########################
    if (mord==1) {
    Bi<-diag(clam0[,1])
    temp4<-cla0^2/di
    As<-diag(temp4[,1])
    mat<-(Wi%*%Bi)-(Wi%*%Mi)%*%As%*%(t(Mi)%*%Wi)
    dinv<-solve(t(z)%*%mat%*%z+U)
    mu0<-exp(x%*%beta_h0 + z%*%v_h0)
    mu<-exp(x%*%beta_h0 + z%*%v_h0)*clam0
    dcla0be<--As%*%(t(Mi)%*%Wi)
    dcla0b<-matrix(0,idx2,p)
    dv_db<-matrix(0,qcum[nrand+1],p)
    xz<-matrix(0,n,p)
    dmu0<-matrix(0,n,p)
    dw_db1<-matrix(0,n,p)
    xk<-matrix(0,n,1)
    ad1<-matrix(0,p,1)
    for (k in 1:p) {
       xk[,1]<-x[,k]
       dv_db[,k] <--dinv%*%(t(z)%*%mat%*%xk)
       xz[,k]<-xk+z%*%(dv_db[,k])
       dcla0b[,k]<-dcla0be%*%(xz[,k])
       dc<-Mi%*%diag(dcla0b[,k])%*%done
       dmu0[,k]<-mu0*(xz[,k])
       dw_db1[,k]<- dc*mu0 + mu*xz[,k]
       temp4<-(2*cla0/di)*(dcla0b[,k])
       dw_db2<-((diag(dmu0[,k]))%*%Mi%*%As%*%t(Mi)%*%Wi)+(Wi%*%Mi%*%diag(temp4[,1])%*%t(Mi)%*%Wi)+(Wi%*%Mi%*%As%*%t(Mi)%*%(diag(dmu0[,k]))) 
       dw_db<-diag(dw_db1[,k])-dw_db2
       ad1[k,1]<-sum(diag(dinv%*%(t(z)%*%dw_db%*%z))) 
    }
    dft1<-dft1-0.5*ad1
   }
#########################################################
    dft2<-t(z)%*%(del-Wi%*%clam0)+(iD%*%oq)-(iD%*%iu_h0) ## gamma frailty
    dft<-rbind(dft1,dft2)
    Bi<-diag(clam0[,1])
    temp4<-cla0^2/di
    Adi<-diag(temp4[,1])
    As<-Adi
    mat<-(Wi%*%Bi)-(Wi%*%Mi)%*%Adi%*%(t(Mi)%*%Wi)
    U<-iD%*%diag(iu_h0[,1]) ## gamma frailty
    H <- rbind(cbind(t(x)%*%mat%*%x, t(x)%*%mat%*%z), cbind(t(z)%*%mat%*%x, t(z)%*%mat%*%z+U))
    Hinv<-solve(H)
    be_h0<- rbind(beta_h0, v_h0)
    be_h <- be_h0 + (Hinv%*%dft)
    beta_h[1:p,1]<-be_h[1:p,1]
    se_beta_h<-matrix(0,p,1)
    for (i in 1:p) se_beta_h[i,1]<-sqrt(Hinv[i,i])
    index2<-qcum[nrand+1]
    index3<-p+1
    index4<-p+qcum[nrand+1]
    v_h[1:index2,1]<-be_h[index3:index4,1]
################################################
    for (i in 1:nrand) {
        ial<-1/alpha_h[i]
        dp<-digamma(ial)
        ddp<-trigamma(ial)
        oq<-matrix(1,q[i],1)
        one<-matrix(1,n,1)
        u_h<-exp(v_h) 
        eta<-x%*%beta_h + z%*%v_h
        expeta<-exp(eta)
        Wi<-diag(expeta[,1])
        Wei<-(Wi%*%Bi)
        U<-iD%*%diag(u_h[,1])
        term<-(t(z)%*%mat%*%z+U)
        C<-matrix(0,qcum[nrand+1],qcum[nrand+1])
        index1<-qcum[i]+1
        index2<-qcum[i+1]
        for (j in index1:index2) C[j,j]<-1
        iA<-iD
        iB<-iA%*%C%*%iA
        c_vh<- iB%*%(u_h-1) ## gamma frailty
        invt<-solve(term)
        dv<-invt%*%c_vh
        dexpeta<-expeta*(z%*%dv)
        dcla0<--(di/((t(Mi)%*%expeta)^2))*(t(Mi)%*%dexpeta)
        temp4<-expeta*(z%*%dv)
        dWi<-diag(temp4[,1])
        dAi<-diag(dcla0[,1])
        temp4<-Mi%*%dAi%*%done
        dBi<-diag(temp4[,1])
        dvec<-2*(cla0*dcla0)
        temp4<-dvec/di
        dAs<-diag(temp4[,1])
        dmat<-(dWi%*%Bi)+(Wi%*%dBi)-(dWi%*%Mi%*%As%*%t(Mi)%*%Wi)-(Wi%*%Mi%*%dAs%*%t(Mi)%*%Wi)-(Wi%*%Mi%*%As%*%t(Mi)%*%dWi)
        uad1<-u_h*dv
        temp4<--iB%*%u_h+ial*uad1
        dia1<-diag(temp4[,1])
        Hd <- rbind(cbind(t(x)%*%dmat%*%x,t(x)%*%dmat%*%z),cbind(t(z)%*%dmat%*%x,t(z)%*%dmat%*%z+ dia1))
    if (dord==0) {
        temp4<--iD%*%iD%*%u_h
        dia1_k<-diag(temp4[,1])
        zero1<-matrix(0,p,p)
        zero2<-matrix(0,p,qcum[nrand+1])
        Hd_k<-rbind(cbind(zero1,zero2),cbind(t(zero2),dia1_k))
        hinv2<-solve(t(z)%*%mat%*%z+U)
        Hd2<-t(z)%*%dmat%*%z+dia1
        dk2<--0.5*sum(diag(Hinv%*%Hd_k))
        vv_h<-matrix(0,q[i],1)
        index3<-1
        index4<-q[i]
        vv_h[1:q[i],1]<-v_h[index1:index2,1]
        uu_h<-exp(vv_h)
        k2<-(t(oq)%*%(vv_h -uu_h))+( q[i]*(-log(alpha_h[i])+1 - dp) )
        mu<-exp(x%*%beta_h)*clam0
        zd<-t(z)%*%del
        zmu<-t(z)%*%mu
        cor1<-(1+(alpha_h[i]*zd))^(-2)
        dk3<-sum(cor1)/12
        k2<--((alpha_h[i]^-2)*k2)+dk2
        dterm<-t(z)%*%dmat%*%z+dia1
        ddv<-- invt%*%(dterm)%*%invt%*%c_vh+invt%*%(-2*ial^3*(u_h-1)+ial^2*uad1 )
        ddcla0<--dAs%*%(t(Mi)%*%Wi%*%z)%*%dv-As%*%(t(Mi)%*%dWi%*%z)%*%dv-As%*%(t(Mi)%*%Wi%*%z)%*%ddv 
        uad2<-(uad1*dv)+(u_h*ddv)
        ddAi<-diag(ddcla0[,1])
        ddclam0<-Mi%*%ddAi%*%done
        ddBi<-diag(ddclam0[,1])
        temp4<-expeta*(z%*%dv)*(z%*%dv)+(expeta*(z%*%ddv))
        ddWi<-diag(temp4[,1])
        ddm1<-(ddWi%*%Bi)+(2*dWi%*%dBi) + (Wi%*%ddBi)
        ddm2<-ddWi%*%Mi%*%As%*%t(Mi)%*%Wi+(2*dWi%*%Mi%*%dAs%*%t(Mi)%*%Wi)+(2*dWi%*%Mi%*%As%*%t(Mi)%*%dWi)
        ddvec<-(2*(dcla0^2)) + (2*(cla0*ddcla0))
        temp4<-ddvec/di
        ddAs=diag(temp4[,1])
        ddm3<-(Wi%*%Mi%*%ddAs%*%t(Mi)%*%Wi)+(2*Wi%*%Mi%*%dAs%*%t(Mi)%*%dWi)+Wi%*%Mi%*%As%*%t(Mi)%*%ddWi
        ddmat<-ddm1-(ddm2+ddm3)
        temp4<-(2*ial^3*u_h)-(2*ial^2*uad1)+(ial*uad2)
        dia2<-diag(temp4[,1])
        Hdd2<-t(z)%*%ddmat%*%z+dia2
        k21<-t(oq)%*%(2*(vv_h -uu_h))+ ( q[i]*( (-2*log(alpha_h[i]))+3 - (2*dp)-((1/alpha_h[i])*ddp))  )
        al<-ial^3
        k21<-al*k21
        cor2<-(1+(alpha_h[i]*zd))^(-3)
        cor22<-cor2*zd
        kcor<-sum(cor22)/6
        kcor<-0
        k22<--k21+0.5*(sum(diag(hinv2*Hdd2))-sum(diag(hinv2*Hd2*hinv2*Hd2)))+kcor
        ialp<-1/k22
        alpha_h[i]<-alpha_h[i] + (ialp*k2)
    }
    if (dord==1 | dord==2) {
        hinv2<-solve(t(z)%*%mat%*%z+U)
        Hd2<-t(z)%*%dmat%*%z+dia1
        dk2<--0.5*sum(diag(Hinv%*%Hd))
        vv_h<-matrix(0,q[i],1)
        index3<-1
        index4<-q[i]
        vv_h[1:index4,1]<-v_h[index1:index2,1]
        uu_h<-exp(vv_h)
        k2<-(t(oq)%*%(vv_h -uu_h))+( q[i]*(-log(alpha_h[i])+1 - dp) )
        mu<-exp(x%*%beta_h)*clam0
        zd<-t(z)%*%del
        zmu<-t(z)%*%mu
        cor1<-(1+(alpha_h[i]*zd))^(-2)
        dk3<-sum(cor1)/12
        if(dord==1) k2<--((alpha_h[i]^(-2))*k2)+dk2
        if(dord==2) k2<--((alpha_h[i]^(-2))*k2)+dk2 +dk3
        dterm<-t(z)%*%dmat%*%z+dia1
        ddv<--invt%*%(dterm)%*%invt%*%c_vh+ invt%*%( -2*ial^3*(u_h-1)+ial^2*uad1 )
        ddcla0<--dAs%*%(t(Mi)%*%Wi%*%z)%*%dv-As%*%(t(Mi)%*%dWi%*%z)%*%dv-As%*%(t(Mi)%*%Wi%*%z)%*%ddv 
        uad2<-(uad1*dv)+(u_h*ddv)
        ddAi<-diag(ddcla0[,1])
        ddclam0<-Mi%*%ddAi%*%done
        ddBi<-diag(ddclam0[,1])
        temp4<-expeta*(z%*%dv)*(z%*%dv)+(expeta*(z%*%ddv))
        ddWi<-diag(temp4[,1])
        ddm1<-(ddWi%*%Bi)+(2*dWi%*%dBi) + (Wi%*%ddBi)
        ddm2<-ddWi%*%Mi%*%As%*%t(Mi)%*%Wi+(2*dWi%*%Mi%*%dAs%*%t(Mi)%*%Wi)+(2*dWi%*%Mi%*%As%*%t(Mi)%*%dWi)
        ddvec<-(2*(dcla0^2)) + (2*(cla0*ddcla0))
        temp4<-ddvec/di
        ddAs<-diag(temp4[,1])
        ddm3<-(Wi%*%Mi%*%ddAs%*%t(Mi)%*%Wi)+(2*Wi%*%Mi%*%dAs%*%t(Mi)%*%dWi)+Wi%*%Mi%*%As%*%t(Mi)%*%ddWi
        ddmat<-ddm1-(ddm2+ddm3)
        temp4<-(2*ial^3*u_h)-(2*ial^2*uad1)+(ial*uad2)
        dia2<-diag(temp4[,1])
        Hdd<-rbind(cbind(t(x)%*%ddmat%*%x,t(x)%*%ddmat%*%z),cbind(t(z)%*%ddmat%*%x,t(z)%*%ddmat%*%z+dia2))
        Hdd2<-t(z)%*%ddmat%*%z+dia2
        k21<-t(oq)%*%(2*(vv_h -uu_h))+(q[i]*((-2*log(alpha_h[i]))+3-(2*dp)-((1/alpha_h[i])*ddp)))
        al<-ial^3
        k21<-al*k21
        cor2<-(1+(alpha_h[i]*zd))^(-3)
        cor22<-cor2*zd
        kcor<-sum(cor22)/6
        if(dord==1) kcor<-0
        k22<--k21+0.5*(sum(diag(Hinv%*%Hdd))-sum(diag(Hinv%*%Hd%*%Hinv%*%Hd)))+kcor
        ialp<-1/k22
        alpha_h[i] <- alpha_h[i] + (ialp*k2)
        if (alpha_h[i]<=0.0) alpha_h[i]<-alpha_h0/2
    }
    }
    res<-list(x,z,y,del,Mi,idx2,t2,di,beta_h0,v_h0,beta_h,v_h,alpha_h0,alpha_h,dft,Hinv,clam0,H,mat,se_beta_h,U)
    return(res)
}

PGFrailty_SE.h<-function(res1,nrand,q,qcum,dord=1) {
x<-res1[1][[1]]
z<-res1[2][[1]]
y<-res1[3][[1]]
del<-res1[4][[1]]
Mi<-res1[5][[1]]
idx2<-res1[6][[1]]
t2<-res1[7][[1]]
di<-res1[8][[1]]
beta_h<-res1[9][[1]]
v_h<-res1[10][[1]]
beta_h1<-res1[11][[1]]
v_h1<-res1[12][[1]]
alpha_h<-res1[13][[1]]
alpha_h1<-res1[14][[1]]
dft<-res1[15][[1]]
Hinv<-res1[16][[1]]
clam0<-res1[17][[1]]
H<-res1[18][[1]]
mat<-res1[19][[1]]

################################################
######## SE for frailty parameter ###############
################################################
    n<-nrow(x)
    p<-ncol(x)
    u_h1<-exp(v_h1)
    mat11<-t(x)%*%mat%*%x
    mat12<-t(x)%*%mat%*%z
    mat13<-t(z)%*%mat%*%z
    oq<-matrix(1,qcum[nrand+1],1)
    oq1<-matrix(1,qcum[nrand+1],1)
    for (i in 1:nrand) {
       index1<-qcum[i]+1
       oq1[index1:qcum[i+1]]<-alpha_h1[i]
    } 
    D<-diag(oq1[,1])
    iD<-solve(D)
    U <- iD%*%diag(u_h1[,1])
    mmat<-mat11-mat12%*%solve(mat13+U)%*%t(mat12)
    hminv<-solve(mmat)
    done<-matrix(1,idx2,1)
    muh<-x%*%beta_h1 + z%*%v_h1
    expeta<-exp(muh)
    Wi<-diag(expeta[,1])
    cla0<-di/(t(Mi)%*%expeta)
    Ai<-diag(cla0[,1])
    temp4<-cla0^2/di
    As<-diag(temp4[,1])
    Bi<-diag(clam0[,1])
    mat<-(Wi%*%Bi)-(Wi%*%Mi)%*%As%*%(t(Mi)%*%Wi)
    Dinv0<-solve(t(z)%*%mat%*%z+U)
    se_lam<-matrix(0,nrand,1)
    for (i in 1:nrand){
    ial<-1/alpha_h1[i]
    index1<-qcum[i]+1
    index2<-qcum[i+1]
    vv_h1<-matrix(0,q[i],1)
    uu_h1<-matrix(0,q[i],1)
    vv_h1[1:q[i],1]<-v_h1[index1:index2,1]
    uu_h1[1:q[i],1]<-u_h1[index1:index2,1]
    c_vh <- ial^2*(uu_h1-1)
    dv<-solve(t(z)%*%mat%*%z+U)%*%c_vh
    dexpeta<-expeta*(z%*%dv)
    dcla0<--(di/((t(Mi)%*%expeta)^2))*(t(Mi)%*%dexpeta)
    dWi<-diag(dexpeta[,1])
    dAi<-diag(dcla0[,1])
    temp4<-Mi%*%dAi%*%done
    dBi<-diag(temp4[,1])
    dvec<-2*(cla0*dcla0)
    temp4<-dvec/di
    dAs<-diag(temp4[,1])
    dmat<-(dWi%*%Bi)+(Wi%*%dBi)-(dWi%*%Mi%*%As%*%t(Mi)%*%Wi)-(Wi%*%Mi%*%dAs%*%t(Mi)%*%Wi)-(Wi%*%Mi%*%As%*%t(Mi)%*%dWi)
    uad1<-uu_h1*dv
    temp4<-(-ial^2*uu_h1)+(ial*uad1) 
    dia_d<-diag(temp4[,1])
    Hd <- rbind(cbind(t(x)%*%dmat%*%x,t(x)%*%dmat%*%z),cbind(t(z)%*%dmat%*%x,t(z)%*%dmat%*%z+ dia_d))
    term<-(t(z)%*%mat%*%z+U)
    invt<-solve(term)
    dterm<-t(z)%*%dmat%*%z + dia_d
    ddv<- - invt%*%(dterm)%*%invt%*%c_vh+ invt%*%( -2*ial^3*(uu_h1-1)+ ial^2*uad1 )
    ddcla0<- -dAs%*%(t(Mi)%*%Wi%*%z)%*%dv-As%*%(t(Mi)%*%dWi%*%z)%*%dv-As%*%(t(Mi)%*%Wi%*%z)%*%ddv 
    ddAi<-diag(ddcla0[,1])
    ddclam0<-Mi%*%ddAi%*%done
    ddBi<-diag(ddclam0[,1])
    temp4<-expeta*(z%*%dv)*(z%*%dv)+(expeta*(z%*%ddv))
    ddWi<-diag(temp4[,1])
    ddm1<-(ddWi%*%Bi)+(2*dWi%*%dBi) + (Wi%*%ddBi)
    ddm2<- ddWi%*%Mi%*%As%*%t(Mi)%*%Wi+(2*dWi%*%Mi%*%dAs%*%t(Mi)%*%Wi)+(2*dWi%*%Mi%*%As%*%t(Mi)%*%dWi)
    ddvec<-(2*(dcla0^2)) + (2*(cla0*ddcla0))
    temp4<-ddvec/di
    ddAs<-diag(temp4[,1])
    ddm3<-(Wi%*%Mi%*%ddAs%*%t(Mi)%*%Wi)+(2*Wi%*%Mi%*%dAs%*%t(Mi)%*%dWi)+Wi%*%Mi%*%As%*%t(Mi)%*%ddWi
    ddmat<-ddm1-(ddm2+ddm3)
    uad2<-(uad1*dv)+(u_h1*ddv)
    temp4<- (2*ial^3*u_h1)-(2*ial^2*uad1)+(ial*uad2) 
    dia_dd<-diag(temp4[,1])
    Hdd <-rbind(cbind(t(x)%*%ddmat%*%x,t(x)%*%ddmat%*%z),cbind(t(z)%*%ddmat%*%x,t(z)%*%ddmat%*%z+ dia_dd))
    H <- rbind(cbind(t(x)%*%mat%*%x,t(x)%*%mat%*%z),cbind(t(z)%*%mat%*%x,t(z)%*%mat%*%z+U))
    Hinv<-solve(H)
    oq<-matrix(1,q[i],1)
    dp<-digamma(ial)
    ddp<-trigamma(ial)
    k21a<-t(oq)%*%(2*(vv_h1 -uu_h1))+( q[i]*( (-2*log(alpha_h1[i]))+3 -(2*dp)-((1/alpha_h1[i])*ddp))  )
    k21a<-(ial^3)*k21a
    d2halp<--k21a- t(oq)%*%(c_vh*dv)
    adalp<-0.5*sum(diag(-Hinv%*%Hd%*%Hinv%*%Hd+ Hinv%*%Hdd))
    dalp_2<-d2halp+adalp
    se_lam[i]<-sqrt(1/dalp_2)
    }
    u_h1<-exp(v_h1)
    U <- iD%*%diag(u_h1[,1])
    oq<-matrix(1,qcum[nrand+1],1)
    one<-matrix(1,n,1)
    zd<-t(z)%*%del
    eta<-x%*%beta_h1 + z%*%v_h1
    expeta<-exp(eta)
    term0<-t(Mi)%*%expeta
    non<-t(one)%*%del
    done<-matrix(1,idx2,1)
    hlike1<-(t(one)%*%(del*eta) )-( t(done)%*%(di*log(term0)) )
    hlike2<-0
    for (i in 1:nrand) {
       oq<-matrix(1,q[i],1)
       index1<-qcum[i]+1
       index2<-qcum[i+1]
       vv_h1<-matrix(0,q[i],1)
       uu_h1<-matrix(0,q[i],1)
       vv_h1[1:q[i],1]<-v_h1[index1:index2,1]
       uu_h1[1:q[i],1]<-u_h1[index1:index2,1]
       i_alp1<-1/alpha_h1[i]
       c_alp1<--log(gamma(i_alp1))-(i_alp1*log(alpha_h1[i]))
       hlike2<- hlike2+t(oq)%*%( (vv_h1-uu_h1)/alpha_h1[i] + c_alp1)
    }
    hlike<-hlike1+hlike2
    pi<-3.14159265359
    H22<-t(z)%*%mat%*%z+U
    zd<-t(z)%*%del
    if (dord==2) {
       temp4<-1/(i_alp1+zd)
       secd<-diag(temp4[,1])
       second<-sum(diag(secd))/12
    } else second<-0
    pvhs<-hlike-0.5*log(det(H22/(2*pi)))
    svhs<-pvhs+second
    adj1<-( (0.5*(p+qcum[nrand+1]))*log(2*pi)) + (0.5*log(det(Hinv)) )
    hpn1<-hlike+ adj1
    hpn2<-pvhs
    hpn3<-hpn1+second
    res<-list(se_lam,hlike,hpn1,hpn2,hpn3)
    return(res)
}


PNFrailty_SE.h<-function(res1,nrand,q,qcum,dord=1) {
x<-res1[1][[1]]
z<-res1[2][[1]]
y<-res1[3][[1]]
del<-res1[4][[1]]
Mi<-res1[5][[1]]
idx2<-res1[6][[1]]
t2<-res1[7][[1]]
di<-res1[8][[1]]
beta_h<-res1[9][[1]]
v_h<-res1[10][[1]]
beta_h1<-res1[11][[1]]
v_h1<-res1[12][[1]]
alpha_h<-res1[13][[1]]
alpha_h1<-res1[14][[1]]
dft<-res1[15][[1]]
Hinv<-res1[16][[1]]
clam0<-res1[17][[1]]
H<-res1[18][[1]]
mat<-res1[19][[1]]
U<-res1[21][[1]]

################################################
######## SE for frailty parameter ###############
################################################
    n<-nrow(x)
    p<-ncol(x)
    u_h1<-exp(v_h1)
    oq<-matrix(1,qcum[nrand+1],1)
    oq1<-matrix(1,qcum[nrand+1],1)
    for (i in 1:nrand) {
       index1<-qcum[i]+1
       oq1[index1:qcum[i+1]]<-alpha_h1[i]
    } 
    D<-diag(oq1[,1])
    iD<-solve(D)
    iA<-iD
    Bi<-diag(clam0[,1])
    muh<-x%*%beta_h1 + z%*%v_h1
    expeta<-exp(muh)
    cla0<-di/(t(Mi)%*%expeta)
    temp4<-cla0^2/di
    As<-diag(temp4[,1])
    Wi<-diag(expeta[,1])
    done<-matrix(1,idx2,1)
    H22<-solve(t(z)%*%mat%*%z+U)
    Hessian<-matrix(0,nrand,nrand)
    for (i in 1:nrand) {
        C<-matrix(0,qcum[nrand+1],qcum[nrand+1])
        index1<-qcum[i]+1
        index2<-qcum[i+1]
        for (j in index1:index2) C[j,j]<-1
        iB1<-iA%*%C%*%iA
        c_vh1<-iB1%*%v_h1
        dv1<-H22%*%c_vh1
        dexpeta1<-expeta*(z%*%dv1)
        dcla01<--(di/((t(Mi)%*%expeta)^2))*(t(Mi)%*%dexpeta1)
        dWi1<-diag(dexpeta1[,1])
        dAi1<-diag(dcla01[,1])
        temp4<-Mi%*%dAi1%*%done
        dBi1<-diag(temp4[,1])
        dvec1<-2*(cla0*dcla01)
        temp4<-dvec1/di
        dAs1<-diag(temp4[,1])
        dmat1<-(dWi1%*%Bi)+(Wi%*%dBi1)-(dWi1%*%Mi%*%As%*%t(Mi)%*%Wi)-(Wi%*%Mi%*%dAs1%*%t(Mi)%*%Wi)-(Wi%*%Mi%*%As%*%t(Mi)%*%dWi1)
        dia1<--iB1
        Hd1 <- rbind(cbind(t(x)%*%dmat1%*%x,t(x)%*%dmat1%*%z),cbind(t(z)%*%dmat1%*%x,t(z)%*%dmat1%*%z+ dia1))
        ddk1<- -0.5*sum(diag(iA%*%C%*%iA%*%C)) +t(v_h1)%*%(iA%*%C%*%iB1)%*%v_h1 -t(dv1)%*%iB1%*%v_h1
        dia11<-(iB1%*%C%*%iA+iA%*%C%*%iB1)
        dv11<--H22%*%((t(z)%*%dmat1%*%z+dia1)%*%dv1-iB1%*%dv1 + dia11%*%v_h1)
        temp4<-(z%*%dv1)*(z%*%dv1)*expeta  +(z%*%dv11)*expeta
        ddW11<-diag(temp4[,1]) 
        ddcla011<- -( dAs1%*%(t(Mi)%*%Wi%*%z)%*%dv1 +As%*%(t(Mi)%*%dWi1%*%z)%*%dv1 +As%*%(t(Mi)%*%Wi%*%z)%*%dv11)  
        temp4<-Mi%*%diag(ddcla011[,1])%*%done
        ddB11<-diag(temp4[,1])
        temp4<-(2*(dcla01^2) + 2*(cla0*ddcla011) ) /di
        ddAs11<-diag(temp4[,1])
        ddm1_11<-(ddW11%*%Bi)+ (2*dWi1%*%dBi1) + (Wi%*%ddB11)
        ddm2_11<- (ddW11%*%Mi%*%As%*%t(Mi)%*%Wi) +(2*dWi1%*%Mi%*%dAs1%*%t(Mi)%*%Wi) +(2*dWi1%*%Mi%*%As%*%t(Mi)%*%dWi1)
        ddm3_11<-(Wi%*%Mi%*%ddAs11%*%t(Mi)%*%Wi) +(2*Wi%*%Mi%*%dAs1%*%t(Mi)%*%dWi1)  +(Wi%*%Mi%*%As%*%t(Mi)%*%ddW11)
        ddmat11<-ddm1_11-(ddm2_11+ddm3_11)
        Hd11 <-rbind(cbind(t(x)%*%ddmat11%*%x,t(x)%*%ddmat11%*%z),cbind(t(z)%*%ddmat11%*%x,t(z)%*%ddmat11%*%z+ dia11 ))
        ddk1<-ddk1 +0.5*sum(diag(-Hinv%*%Hd1%*%Hinv%*%Hd1+ Hinv%*%Hd11))
        Hessian[i,i]<-ddk1
        for (kk in 1:nrand) {
          if (kk>i) {
             D<-matrix(0,qcum[nrand+1],qcum[nrand+1])
             index1<-qcum[kk]+1
             index2<-qcum[kk+1]
             for (j in index1:index2) D[j,j]<-1
             iB2<-iA%*%D%*%iA
             c_vh2<-iB2%*%v_h1
             dv2<-H22%*%c_vh2
             dexpeta2<-expeta*(z%*%dv2)
             dcla02<--(di/((t(Mi)%*%expeta)^2))*(t(Mi)%*%dexpeta2)
             dWi2<-diag(dexpeta2[,1])
             dAi2<-diag(dcla02[,1])
             temp4<-Mi%*%dAi2%*%done
             dBi2<-diag(temp4[,1])
             dvec2<-2*(cla0*dcla02)
             temp4<-dvec2/di
             dAs2<-diag(temp4[,1])
             dd12<--0.5*sum(diag(iA%*%D%*%iA%*%C))+0.5*t(v_h1)%*%(iA%*%D%*%iB1+iA%*%C%*%iB2)%*%v_h1-t(dv1)%*%iB2%*%v_h1
             dia12<-(iB1%*%D%*%iA+iA%*%D%*%iB1)
             dv12<--H22%*%((t(z)%*%dmat1%*%z+dia1)%*%dv2-iB2%*%dv1 + dia12%*%v_h1)
             temp4<-(z%*%dv1)*(z%*%dv2)*expeta + (z%*%dv12)*expeta
             ddW12<-diag(temp4[,1])
             ddcla012<--( dAs2%*%(t(Mi)%*%Wi%*%z)%*%dv1 +As%*%(t(Mi)%*%dWi2%*%z)%*%dv1 +As%*%(t(Mi)%*%Wi%*%z)%*%dv12)  
             temp4<-Mi%*%diag(ddcla012[,1])%*%done
             ddB12<-diag(temp4[,1])
             temp4<-(2*(dcla02*dcla01) + 2*(cla0*ddcla012) )/di
             ddAs12<-diag(temp4[,1])
             ddm1_12<-(ddW12%*%Bi)+ (dWi1%*%dBi2 + dWi2%*%dBi1) + (Wi%*%ddB12)
             ddm2_12<-(ddW12%*%Mi%*%As%*%t(Mi)%*%Wi) +(dWi1%*%Mi%*%dAs2%*%t(Mi)%*%Wi) +(dWi1%*%Mi%*%As%*%t(Mi)%*%dWi2)
             ddm3_12<-(dWi2%*%Mi%*%dAs1%*%t(Mi)%*%Wi) +(Wi%*%Mi%*%ddAs12%*%t(Mi)%*%Wi)  +(Wi%*%Mi%*%dAs1%*%t(Mi)%*%dWi2)
             ddm4_12<-(dWi2%*%Mi%*%As%*%t(Mi)%*%dWi1) +(Wi%*%Mi%*%dAs2%*%t(Mi)%*%dWi1)  +(Wi%*%Mi%*%As%*%t(Mi)%*%ddW12)
             ddmat12<-ddm1_12-(ddm2_12+ddm3_12+ddm4_12)
             dmat2<-(dWi2%*%Bi)+(Wi%*%dBi2)-(dWi2%*%Mi%*%As%*%t(Mi)%*%Wi)-(Wi%*%Mi%*%dAs2%*%t(Mi)%*%Wi)-(Wi%*%Mi%*%As%*%t(Mi)%*%dWi2)
             dia2<--iB2
             Hd2 <- rbind(cbind(t(x)%*%dmat2%*%x,t(x)%*%dmat2%*%z),cbind(t(z)%*%dmat2%*%x,t(z)%*%dmat2%*%z+ dia2))
             Hd12<-rbind(cbind(t(x)%*%ddmat12%*%x,t(x)%*%ddmat12%*%z),cbind(t(z)%*%ddmat12%*%x,t(z)%*%ddmat12%*%z+ dia12))
             dd12<-dd12 +0.5*sum(diag(-Hinv%*%Hd1%*%Hinv%*%Hd2+ Hinv%*%Hd12))
             Hessian[i,kk]<-Hessian[kk,i]<-dd12
         }
       }
     }
     iAp<-solve(Hessian)
     se_lam<-sqrt(diag(iAp))
     eta<-x%*%beta_h1 + z%*%v_h1
     expeta<-exp(eta)
     one<-matrix(1,n,1)
     done<-matrix(1,idx2,1)
     oq<-matrix(1,qcum[nrand+1],1)
     pi<-3.14159265359
     term0<-t(Mi)%*%expeta
     hlike1<-(t(one)%*%(del*eta) )-( t(done)%*%(di*log(term0)))
     hlike2<-0
     hlike3<-0
     for (i in 1:nrand) {
         hlike2<-hlike2-(q[i]/2)*log(2*pi)-( (q[i]/2)*log(alpha_h1[i]))
         index1<-qcum[i]+1
         index2<-qcum[i+1]
         vv_h1<-matrix(0,q[i],1)
         vv_h1[1:q[i],1]<-v_h1[index1:index2,1]
         hlike3<-hlike3-(t(vv_h1)%*%vv_h1)/(2*alpha_h1[i])
     }
     hliken<-hlike1+hlike2+hlike3
     adj1<-( (0.5*(p+qcum[nrand+1]))*log(2*pi)) + (0.5*log(det(Hinv)) )
     hpn1<-hliken+ adj1
     muu<-exp(x%*%beta_h1)*clam0
     zmu<-t(z)%*%muu
     u_h1<-exp(v_h1)
     second<-0
     for (i in 1:nrand) {
         ialph1<-1/alpha_h1[i]
         a21<-(zmu*u_h1)+ialph1
         b31<-zmu*u_h1
         S11<-3*(b31/(a21^2))
         S21<-5*((b31^2)/(a21^3))
         temp4<-S11-S21
         S31<-diag(temp4[,1])
         second<-second-sum(diag(S31))/24
    }
    H22<-t(z)%*%mat%*%z+iD
    hpn2<-hliken-0.5*log(det(H22/(2*pi)))
    hpn3<-hpn1+second
    res<-list(se_lam,hliken,hpn1,hpn2,hpn3)
    return(res)
}


subbars<-function (term) 
{
    if (is.name(term) || !is.language(term)) 
        return(term)
    if (length(term) == 2) {
        term[[2]] <- subbars(term[[2]])
        return(term)
    }
    stopifnot(length(term) >= 3)
    if (is.call(term) && term[[1]] == as.name("|")) 
        term[[1]] <- as.name("+")
    for (j in 2:length(term)) term[[j]] <- subbars(term[[j]])
    term
}

nobars<-function (term) 
{
    if (!("|" %in% all.names(term))) 
        return(term)
    if (is.call(term) && term[[1]] == as.name("|")) 
        return(NULL)
    if (length(term) == 2) {
        nb <- nobars(term[[2]])
        if (is.null(nb)) 
            return(NULL)
        term[[2]] <- nb
        return(term)
    }
    nb2 <- nobars(term[[2]])
    nb3 <- nobars(term[[3]])
    if (is.null(nb2)) 
        return(nb3)
    if (is.null(nb3)) 
        return(nb2)
    term[[2]] <- nb2
    term[[3]] <- nb3
    term
}

expandSlash<-function (bb) 
{
    if (!is.list(bb)) 
        return(expandSlash(list(bb)))
    unlist(lapply(bb, function(x) {
        if (length(x) > 2 && is.list(trms <- slashTerms(x[[3]]))) 
            return(lapply(unlist(makeInteraction(trms)), function(trm) substitute(foo | 
                bar, list(foo = x[[2]], bar = trm))))
        x
    }))
}

findbars<-function(term) 
{
    if (is.name(term) || !is.language(term)) 
        return(NULL)
    if (term[[1]] == as.name("(")) 
        return(findbars(term[[2]]))
    if (!is.call(term)) 
        stop("term must be of class call")
    if (term[[1]] == as.name("|")) 
        return(term)
    if (length(term) == 2) 
        return(findbars(term[[2]]))
    c(findbars(term[[2]]), findbars(term[[3]]))
}

##### Check if there were two terms entered in a correlated fashion #####
##### Current software does not allow for correlated random effects #####
findplus<-function(term)
{
    if (is.numeric(term))
        return(0)        
    if (!is.language(term)) 
        return(NULL)
    if (length(term)==1) return(0)
    if (term[[1]] == as.name("|"))
        return(findplus(term[[2]]))
    if (!is.call(term))
        stop("term must be of class call")
    if (term[[1]] == as.name("+"))
        return(1)
    if (term[[1]] == as.name("-"))
        return(-1)
}
    
slashTerms<-function (x) 
{
    if (!("/" %in% all.names(x))) 
        return(x)
    if (x[[1]] != as.name("/")) 
        stop("unparseable formula for grouping factor")
    list(slashTerms(x[[2]]), slashTerms(x[[3]]))
}

makeInteraction <- function(x)
### from a list of length 2 return recursive interaction terms
{
    if (length(x) < 2) return(x)
    trm1 <- makeInteraction(x[[1]])
    trm11 <- if(is.list(trm1)) trm1[[1]] else trm1
    list(substitute(foo:bar, list(foo=x[[2]], bar = trm11)), trm1)
}
