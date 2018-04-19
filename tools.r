is_true = function(x) {
    if (is.character(x)) {
        return(grepl("<", x))
    }
    else {
        return(x<=0.05)
    }
}

judge_pval = function(a,b,c){
    out = 0
    for (aa in list(a,b,c)){
        out = out + is_true(aa) 
    }
    return(out)
}

set_ranges = function( gr, nm='' ){
   
    upgr=flank(gr, width=100000, start=T)
    dngr=flank(gr, width=100000, start=F)
    newgr=flank(gr, width=100000, both=T)
    end(newgr) <- end(newgr)+width(gr)
    grset = list(gr, trim(upgr), trim(dngr), trim(newgr))
    names(grset) = paste(c("x", "up", "dn", "all"), rep(nm, 4), sep="")
    return(grset)
    
}

prepare_analysis_x = function( xgr ){
    xset1 = set_ranges(xgr) 
    xset2 = set_ranges(xgr[grepl("Nrps", xgr$values),], "nrps")
    xset3 = set_ranges(xgr[grepl("T1pks", xgr$values),], "t1pks")
    return(c(xset1, xset2, xset3))
}

prepare_analysis_y = function( ygr ){
    liney=ygr[ygr$values=="LINE",]
    ltry=ygr[ygr$values=="LTR",]
    dnay=ygr[ygr$values=="DNA",]
    lly=ygr[ygr$values=="LINE" | ygr$values=="LTR", ]
    ally = list(ygr, liney, ltry, dnay, lly)
    names(ally) = c("y", "line", "ltr", "dna", "ll")
    return(ally)
}

get_row_value = function( xobj ) {
    routp = c(0,0,0,0)
    routp[1]= xobj$query.population
    routp[2]= xobj$reference.population
    routp[3]= sum(xobj$projection.test.lower.tail, xobj$scaled.absolute.min.distance.sum.lower.tail, xobj$jaccard.measure.lower.tail)
    routp[4]= judge_pval(xobj$projection.test.p.value, xobj$jaccard.measure.p.value, xobj$scaled.absolute.min.distance.sum.p.value)
    return( routp )
}

single_analysis = function( xx, yy ){
    y_to_x = GenometriCorrelation(query=yy, reference=xx)
    ss = length(y_to_x)
    output = matrix(0L, nrow=ss, ncol=4)
    rownames(output)=names(y_to_x)
    colnames(output)=c("Repeats", "GC", "AveLowTail", "AvePval")
    for (mm in 1:ss ) {
        output[mm,] = get_row_value( y_to_x[[mm]] )
        }
    return(output)
}

geometricorr_analysis = function(allx, ally) {
    namesx = names(allx)
    namesy = names(ally)
    res = c()
    idx = 1
    for (nmx in namesx) {
        for (nmy in namesy) {
            xx = allx[[nmx]]
            yy = ally[[nmy]]
            output <- single_analysis(xx,yy)
            nm = paste(nmx, nmy, sep="-")
            res[[idx]] = output
            names(res)[idx] = nm
            idx = idx+1
        }
    }
    return( res )
}

analysis_gc = function( xgr, ygr ){

    allx = prepare_analysis_x( xgr )
    #names(allx) = c("x", "up", "dn", "uptodn","nrps", "nrpsup", "nrpsdn", "nrpsall", "t1pks", "t1pksup","t1pksdn","t1pksall")
    ally = prepare_analysis_y( ygr )
    
    res = geometricorr_analysis( allx, ally )
    
    return(res)
}
