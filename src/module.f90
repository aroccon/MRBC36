module velocity
    double precision, allocatable :: u(:,:),  v(:,:), rhsu(:,:), rhsv(:,:), rhsu_o(:,:), rhsv_o(:,:)
    double precision, allocatable :: rhsp(:,:),  p(:,:), pext(:,:)
    double complex, allocatable :: rhspc(:,:), pc(:,:), rhs(:)
end module velocity

module phase
    double precision, allocatable :: rhsphi(:,:), phi(:,:), normx(:,:), normy(:,:), fxst(:,:), fyst(:,:)
end module phase

module temperature
    double precision, allocatable :: rhstemp(:,:), temp(:,:)
end module temperature