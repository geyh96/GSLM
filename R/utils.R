#' The function postpropose the interaction selection problem
#' 
#' @param bool_mat_with_name translate the bool matrix to selected interaction term

#' @return return a vector  contains the interaction selected.
#' @export
get_order2_vec=function(bool_mat_with_name){
#  bool_mat_with_name= result_IPDC>Kappa_IPDC$Thres_Selected
  # bool_mat_with_name=DSL_result1_interaction>Kappa_result1$Thres_Selected
        mat_Ind = which(bool_mat_with_name==1 ,   arr.ind=T)
        rowname=rownames(bool_mat_with_name)
        colname=colnames(bool_mat_with_name)
        n_selected = dim(mat_Ind)[1]
        temp=c()
        
        if(n_selected>0){
          for(ii in 1:n_selected){
            temp=c(temp,paste_with_order(rowname[mat_Ind[ii,1]],colname[mat_Ind[ii,2]]))
          }
        return(unique(temp))
        }
        if(n_selected==0){
          # temp=c()
          return(temp)
        }      
}

#' The function postpropose the interaction selection problem
#' 
#' @param str a str like "X1X2"
#' @return return a vector  contains the interaction selected.
str_to_var=function(str){
    # str="X11"
    return(as.numeric(strsplit(str,split="X",fixed=T)[[1]][2]))
  }
#' The function postpropose the interaction selection problem
#' 
#' @param str1 a str like "X1"
#' @param str2 a str like "X2"
#' @return return a vector  like "X1X2"
paste_with_order=function(str1,str2){
    num1 = str_to_var(str1)
    num2 = str_to_var(str2)
    if(num1<num2){
      result = paste(str1,str2,sep="")
    }
    if(num1>num2){
      result = paste(str2,str1,sep="")
    }
    return(result)
}