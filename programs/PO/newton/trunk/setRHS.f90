subroutine setRHS(yi,y,RHS)

real(dp),  dimension(:), intent(in) :: yi, y
real(dp),  dimension(:,:), intent(out) :: RHS

RHS(1:size(RHS)-1) = (y-yi)
RHS(size(RHS)) = 0.0_dp 

end subroutine