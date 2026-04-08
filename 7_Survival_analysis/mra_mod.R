# Versión modificada de corto::mra que acepta un vector como firma de expresión diferencial
mra_mod <- function(expmat1, expmat2 = NULL, regulon, minsize = 10, nperm = NULL, 
                    nthreads = 2, verbose = FALSE, atacseq = NULL) {
  
  if (is.null(nperm)) {
    if (is.null(expmat2)) {
      nperm <- 10
    } else {
      nperm <- 1000
    }
  }
  
  # Filtrar regulones por tamaño mínimo
  regsizes <- sapply(regulon, function(x) length(x$likelihood))
  regulon <- regulon[regsizes >= minsize]
  centroids <- sort(names(regulon))
  targets <- sort(unique(unlist(sapply(regulon, function(x) names(x$tfmode)))))
  
  # Construir matriz de red netmat
  netmat <- matrix(0, nrow = length(centroids), ncol = length(targets),
                   dimnames = list(centroids, targets))
  for (centroid in centroids) {
    vec <- sign(regulon[[centroid]]$tfmode) * regulon[[centroid]]$likelihood
    netmat[centroid, names(vec)] <- vec
  }
  
  # ============================================================
  # NUEVO: Si expmat1 es un vector y expmat2 es NULL -> modo firma
  # ============================================================
  if (is.vector(expmat1) && is.null(expmat2)) {
    if (verbose) message("Modo firma: usando vector proporcionado como expresión diferencial")
    
    sig <- expmat1
    # Filtrar genes que estén en la red
    common <- intersect(names(sig), colnames(netmat))
    if (length(common) == 0) stop("No hay genes comunes entre la firma y la red")
    if (verbose) message(length(common), " genes comunes entre firma y red")
    sig <- sig[common]
    netmat <- netmat[, common, drop = FALSE]
    
    # Normalizar netmat (igual que original)
    netmat <- apply(netmat, 2, function(x) x / (sum(x != 0)^0.5))
    
    # Calcular scores observados
    scores <- (netmat %*% sig)[, 1]
    
    # Función de permutación nula: permutar los nombres del vector sig
    nullsig_perm <- function(seed) {
      set.seed(seed)
      nullsig <- sample(sig, length(sig), replace = FALSE)
      names(nullsig) <- names(sig)
      nullscores <- (netmat %*% nullsig)[, 1]
      return(nullscores)
    }
    
    # Generar permutaciones en paralelo
    cl <- parallel::makeCluster(nthreads)
    nullscores <- pbapply::pbsapply(cl = cl, X = 1:nperm, FUN = nullsig_perm)
    parallel::stopCluster(cl)
    
    # Calcular NES y p-valores (igual que original)
    nes <- apply(cbind(scores, nullscores), 1, function(x) {
      myscore <- x[1]
      morescores <- x[2:length(x)]
      mu <- mean(morescores)
      sigma <- sd(morescores)
      p <- pnorm(abs(myscore), mean = mu, sd = sigma, lower.tail = FALSE) * 2
      if (p == 0) p <- .Machine$double.xmin
      mynes <- p2z(p) * sign(myscore)
      if (myscore == 0) mynes <- 0
      return(mynes)
    })
    
    outlist <- list(nes = nes, pvalue = z2p(nes), sig = sig, regulon = regulon)
    return(outlist)
  }
  
  # ============================================================
  # Comportamiento original con dos matrices (expmat2 no NULL)
  # ============================================================
  if (!is.null(expmat2)) {
    # ... (todo el código original que ya tenías para dos matrices) ...
    # Lo copio íntegro para que la función funcione también con dos matrices
    vargenes <- apply(expmat1, 1, var)
    expmat11 <- expmat1[which(vargenes > 0), ]
    rm(vargenes)
    vargenes <- apply(expmat2, 1, var)
    expmat22 <- expmat2[which(vargenes > 0), ]
    rm(vargenes)
    if (nrow(expmat11) < nrow(expmat1)) {
      message("Removed ", nrow(expmat1) - nrow(expmat11), " rows with zero variance")
    }
    if (nrow(expmat22) < nrow(expmat2)) {
      message("Removed ", nrow(expmat2) - nrow(expmat22), " rows with zero variance")
    }
    common <- intersect(rownames(expmat11), rownames(expmat22))
    expmat1 <- expmat11[common, ]
    expmat2 <- expmat22[common, ]
    rm(expmat11, expmat22)
    sig <- setNames(apply(cbind(expmat1, expmat2), 1, function(x) {
      x1 <- x[1:ncol(expmat1)]
      x2 <- x[(ncol(expmat1) + 1):length(x)]
      tt <- t.test(x1, x2)
      return(tt$statistic)
    }), rownames(expmat1))
    sig <- sig[!is.na(sig)]
    common <- intersect(names(sig), colnames(netmat))
    sig <- sig[common]
    netmat <- netmat[, common]
    netmat <- apply(netmat, 2, function(x) x / (sum(x != 0)^0.5))
    scores <- (netmat %*% sig)[, 1]
    
    nullsigperm1 <- function(seed = 0, expmat1, expmat2, netmat) {
      permat <- cbind(expmat1, expmat2)
      set.seed(seed)
      permat <- permat[sample(nrow(permat)), ]
      permat <- permat[, sample(ncol(permat))]
      nullsig <- setNames(apply(permat, 1, function(x) {
        x1 <- x[1:ncol(expmat1)]
        x2 <- x[(ncol(expmat1) + 1):length(x)]
        tt <- t.test(x1, x2)
        return(tt$statistic)
      }), rownames(permat))
      nullsig <- nullsig[colnames(netmat)]
      nullscores <- (netmat %*% nullsig)[, 1]
      return(nullscores)
    }
    
    cl <- parallel::makeCluster(nthreads)
    nullscores <- pbapply::pbsapply(cl = cl, X = 1:nperm, FUN = nullsigperm1,
                                    expmat1 = expmat1, expmat2 = expmat2, netmat = netmat)
    parallel::stopCluster(cl)
    
    nes <- apply(cbind(scores, nullscores), 1, function(x) {
      myscore <- x[1]
      morescores <- x[2:length(x)]
      mu <- mean(morescores)
      sigma <- sd(morescores)
      p <- pnorm(abs(myscore), mean = mu, sd = sigma, lower.tail = FALSE) * 2
      if (p == 0) p <- .Machine$double.xmin
      mynes <- p2z(p) * sign(myscore)
      if (myscore == 0) mynes <- 0
      return(mynes)
    })
    
    outlist <- list(nes = nes, pvalue = z2p(nes), sig = sig, regulon = regulon)
    return(outlist)
  }
  
  # Si no se cumplió ninguna condición (expmat1 no es vector y expmat2 es NULL) -> error
  stop("Debe proporcionar dos matrices (expmat1 y expmat2) o un vector como firma (expmat1 vector, expmat2 NULL)")
}

# Funciones auxiliares (si no existen en el entorno)
if (!exists("p2z")) p2z <- function(p) qnorm(p/2, lower.tail = FALSE)
if (!exists("z2p")) z2p <- function(z) 2 * pnorm(-abs(z))