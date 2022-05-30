''' function to calculate adjusted R2 '''
def f_adjustedR2(Y, X, R2):
    ''' n = number of samples; p = number of predictors; R2 = R^2 '''
    n = len(Y)
    try :
        p = X1r.shape[1]
    except:
        p = 1
    return 1 - (1-R2)*(n-1)/(n-p-1)