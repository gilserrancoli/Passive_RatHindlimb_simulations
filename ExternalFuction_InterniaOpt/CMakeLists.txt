cmake_minimum_required (VERSION 2.8)
set (CMAKE_CXX_STANDARD 11)

set(TEST_TARGET RightRatHindlimb_Zhong_InertiaVar)
# set(Name1 robotarm6dof_3_1_2)
# set(Name2 robotarm6dof_3_1_2_exactHessian)
project(${TEST_TARGET})


add_library(RightRatHindlimb_Zhong_InertiaVar SHARED foo_jac.c)

install(TARGETS RightRatHindlimb_Zhong_InertiaVar 
	PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE 
				GROUP_READ GROUP_WRITE GROUP_EXECUTE 
				WORLD_READ WORLD_EXECUTE
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	)
		