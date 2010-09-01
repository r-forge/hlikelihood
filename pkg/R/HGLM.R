HGLM <-
function(y = NULL, X = NULL, Z = NULL, X.disp = NULL, family = gaussian(link = identity),
random.family = gaussian(link = identity), method = 'EQL', conv = 1e-4, maxit = 20, 
fixed = NULL, random = NULL, disp = NULL, link.disp = 'log', disp.random = NULL,
data = NULL, data.random = NULL, fix.disp = NULL, Offset = NULL, Weights = NULL,
disp.start = 0, binomial.N = NULL, start.fixed = NULL, start.random = NULL, start.disp = NULL,
start.disp.random = NULL, info = TRUE, debug = FALSE, contrasts = NULL) {

	DistResp <- list()
	DistRand <- list()
	Link <- list()
	LapFix <- list()
	ODEst <- list()
	ODEstVal <- list()
	formulaMain <- list()
	formulaOD <- list()
	formulaRand <- list()
	DataMain <- list()
	DataRand <- list()
	BinomialDen <- list()
	StartBeta <- list()
	StartVs <- list()
	StartRGamma <- list()
	INFO <- list()
	DEBUG <- list()
	CONV <- list() 

	# mc <<- match.call()

	if (!sum(method == c('EQL', 'HL01', 'HL11'))) stop('Invalid estimation method!')
	if (method == 'EQL') {
		startval <- c(start.fixed, start.random, start.disp.random, start.disp)
		if (!is.null(Z)) {
			if (!is.null(Offset)) {
				if (!is.null(Weights)) {
					HGLM.val <- hglm(y = y, X = X, Z = Z, X.disp = X.disp, family = family,
								rand.family = random.family, method = method, conv = conv, maxit = maxit, 
								startval = startval, fixed = fixed, random = random, disp = disp, link.disp = link.disp,
								data = data, fix.disp = fix.disp, offset = Offset, weights = Weights)
				} else {
					HGLM.val <- hglm(y = y, X = X, Z = Z, X.disp = X.disp, family = family,
								rand.family = random.family, method = method, conv = conv, maxit = maxit, 
								startval = startval, fixed = fixed, random = random, disp = disp, link.disp = link.disp,
								data = data, fix.disp = fix.disp, offset = Offset)
				}
			} else {
				if (!is.null(Weights)) {
					HGLM.val <- hglm(y = y, X = X, Z = Z, X.disp = X.disp, family = family,
								rand.family = random.family, method = method, conv = conv, maxit = maxit, 
								startval = startval, fixed = fixed, random = random, disp = disp, link.disp = link.disp,
								data = data, fix.disp = fix.disp, weights = Weights)
				} else {
					HGLM.val <- hglm(y = y, X = X, Z = Z, X.disp = X.disp, family = family,
								rand.family = random.family, method = method, conv = conv, maxit = maxit, 
								startval = startval, fixed = fixed, random = random, disp = disp, link.disp = link.disp,
								data = data, fix.disp = fix.disp)
				}
			}
		} else {
			#print('02')
			#print(Weights)
			if (!is.null(Offset)) {
				if (!is.null(Weights)) {
					HGLM.val <- hglm(fixed = fixed, random = random, family = family,
								disp = disp, link.disp = link.disp, data = data, 
								rand.family = random.family, method = method, conv = conv, maxit = maxit,
								startval = startval, fix.disp = fix.disp, offset = Offset, weights = Weights)
				} else {
					HGLM.val <- hglm(fixed = fixed, random = random, family = family,
								disp = disp, link.disp = link.disp, data = data, 
								rand.family = random.family, method = method, conv = conv, maxit = maxit,
								startval = startval, fix.disp = fix.disp, offset = Offset)
				}
			} else {
				if (!is.null(Weights)) {
					HGLM.val <- hglm(fixed = fixed, random = random, family = family,
							disp = disp, link.disp = link.disp, data = data, 
							rand.family = random.family, method = method, conv = conv, maxit = maxit,
							startval = startval, fix.disp = fix.disp, weights = Weights)
				} else {
					HGLM.val <- hglm(fixed = fixed, random = random, family = family,
							disp = disp, link.disp = link.disp, data = data, 
							rand.family = random.family, method = method, conv = conv, maxit = maxit,
							startval = startval, fix.disp = fix.disp)
				}
			}
		}
	} else {
		#print('1')
		if (method == 'HL01') lap.fix <- FALSE else lap.fix <- TRUE
		if (is.character(family)) family <- get(family)
		if (is.function(family)) family <- eval(family)
		if (is.character(random.family)) random.family <- get(random.family)
		if (is.function(random.family)) random.family <- eval(random.family)
		#print('11')
		if (family$family == 'gaussian') DistResp <<- 'Normal'
		if (family$family == 'binomial') DistResp <<- 'Binomial'
		if (family$family == 'poisson') DistResp <<- 'Poisson'
		if (family$family == 'Gamma') DistResp <<- 'Gamma'
		if (random.family$family == 'gaussian') DistRand <<- 'Normal'
		if (random.family$family == 'Beta') DistRand <<- 'Beta'
		if (random.family$family == 'Gamma') DistRand <<- 'Gamma'
		if (random.family$family == 'inverse.gamma') DistRand <<- 'IGamma'
		if (family$link == 'identity') Link <<- 'Identity'
		if (family$link == 'logit') Link <<- 'Logit'
		if (family$link == 'log') Link <<- 'Log'
		if (family$link == 'inverse') Link <<- 'Inverse'
		#print('12')
		LapFix <<- lap.fix
		if (is.null(fix.disp)) ODEst <<- TRUE else ODEst <<- FALSE
		ODEstVal <<- disp.start
		fixed.terms <- row.names(attr(terms(fixed), 'factors'))
		random.terms <- row.names(attr(terms(random), 'factors'))
		fixed.part <- fixed.terms[2]
		random.part <- paste('(', random.terms[1], ')', sep = '')
		if (length(fixed.terms) > 2) for (i in 3:length(fixed.terms)) fixed.part <- paste(fixed.part, fixed.terms[i], sep = ' + ')
		if (length(random.terms) > 1) for (i in 2:length(random.terms)) random.part <- paste(random.part, paste('(', random.terms[i], ')', sep = ''), sep = ' + ')
		formulaMain <<- as.formula(paste(fixed.terms[1], ' ~ ', fixed.part, ' + ', random.part, sep = ''))
		formulaOD <<- disp
		formulaRand <<- disp.random
		DataMain <<- data
		DataRand <<- data.random
		Offset <<- Offset
		BinomialDen <<- binomial.N
		StartBeta <<- start.fixed
		StartVs <<- start.random
		StartRGamma <<- start.disp.random
		INFO <<- info
		DEBUG <<- debug
		#na.action <- na.action
		contrasts <- contrasts
		CONV <<- conv
		#print('2')
		#HGLMMM.env <- environment(HGLMfit)
		# replace HGLMfit default
		n.data.random <- length(data.random)
		DataRandlist <- 'list(DataRand[[1]]'
		if (n.data.random > 1) {
			for (q in 2:n.data.random) {
				DataRandlist <- paste(DataRandlist, ', DataRand[[', q, ']]', sep = '')
			}
		}
		DataRandlist <- paste(DataRandlist, ')', sep = '')
		towrite <- paste('HGLM.val <<- HGLMfit(DistResp = DistResp, DistRand = DistRand, Link = Link, LapFix = LapFix, ODEst = ODEst, ODEstVal = ODEstVal, formulaMain = formulaMain, formulaOD = formulaOD, formulaRand = formulaRand, DataMain = DataMain, DataRand =', DataRandlist, ', Offset = Offset, BinomialDen = BinomialDen, StartBeta = StartBeta, StartVs = StartVs, StartRGamma = StartRGamma, INFO = INFO, DEBUG = DEBUG, contrasts = contrasts, CONV = CONV)')
		write.table(towrite, 'temp.R', col.names = FALSE, row.names = FALSE, quote = FALSE)
    	#appendHGLMfit()
    	#HGLMfit(DistResp = DistResp, DistRand = DistRand, Link = Link, LapFix = LapFix, ODEst = ODEst, ODEstVal = ODEstVal, formulaMain = formulaMain, formulaOD = formulaOD, formulaRand = formulaRand, DataMain = DataMain, DataRand = list(data.random[[1]]) , Offset = Offset, BinomialDen = BinomialDen, StartBeta = StartBeta, StartVs = StartVs, StartRGamma = StartRGamma, INFO = INFO, DEBUG = DEBUG, contrasts = contrasts, CONV = CONV)
    	#print('3')
    	source('temp.R')
    	file.remove('temp.R')
	}
	return(HGLM.val)
}

