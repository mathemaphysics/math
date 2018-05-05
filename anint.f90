subroutine canint(x,y)
    implicit none
    real*8 :: x,y
    y = dnint(x)
end subroutine

subroutine canintf(x,y)
    implicit none
    real*4 :: x,y
    y = anint(x)
end subroutine

