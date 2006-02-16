module linearAlgebra
  
  use definitions
  implicit none 

CONTAINS

  function cross_product(vec1,vec2)

    real(kind=dp), dimension(3) :: cross_product
    real(kind=dp), dimension(3),intent(in) :: vec1, vec2

    cross_product(1) = vec1(2) * vec2(3) - vec1(3) * vec2(2)
    cross_product(2) = vec1(3) * vec2(1) - vec1(1) * vec2(3)
    cross_product(3) = vec1(1) * vec2(2) - vec1(2) * vec2(1)

  end function cross_product

  function length(vec1)

    real(kind=dp) :: length
    real(kind=dp), dimension(3),intent(in) :: vec1

    length = sqrt(dot_product(vec1, vec1))
  end function length
  
end module linearAlgebra
