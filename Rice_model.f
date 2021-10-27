      PROGRAM Rice
      IMPLICIT NONE
      character(len = 40) :: treatment
      real, dimension (125) :: size_small_spot, size_large_spot
      real, dimension (10) :: small_flounder_case_2 ! for case 2 (small flounder only)
      real, dimension (4) :: large_flounder_case_3  ! for case 3 (large flounder only)
      real, dimension (5) :: small_flounder_case_4  ! for case 4 (mixed)
      real, dimension (2) :: large_flounder_case_4  ! for case 4 (mixed)
      integer :: simulation, days, i, j, k, n_spot
      real, parameter :: prob_predation = 0.5
      real, parameter :: ratio_min = 1.9
      real, parameter :: daily_mort = 0.018
      real, parameter :: prob_of_enc_small = 0.1225
      real, parameter :: prob_of_enc_large = 0.2625
      real, parameter :: growth_small_flounder = 0.77 ! mm/day
      real, parameter :: growth_large_flounder = 0.73 ! mm/day
      integer :: num_of_small_flounder
      integer :: num_of_large_flounder
      real(kind = 8) :: size_of_spot
      real(kind = 8) :: size_of_flounder
      real(kind = 8) :: lower_daily_mort
      real(kind = 8) :: upper_daily_mort
      ! declaring the variables for assigning the size distributions to the spot and flounder
      real(kind = 8) :: mu, a, b, sigma
      real(kind = 8) :: mu_small_spot,lower_small_spot,upper_small_spot
      real(kind = 8) :: mu_large_spot,lower_large_spot,upper_large_spot
      real(kind = 8) :: mu_small_flounder
      real(kind = 8) :: lower_small_flounder,upper_small_flounder
      real(kind = 8) :: mu_large_flounder
      real(kind = 8) :: lower_large_flounder,upper_large_flounder
      integer :: seed, num_of_alive, daily_mortality, alive
      real(kind = 8) :: sigma_spot, sigma_flounder
      real(kind = 8) :: min, max
      
      ! Naming the array that will change as spot are removed due to predation and daily mortality
      real(kind = 8), dimension(5,125) :: final_spot
      real(kind = 8), dimension(9,125) :: small_spot, large_spot
      real(kind = 8), dimension(5,10) :: small_flounder, large_flounder

      ! declaring the variables for the probability of an encounter and predation occuring
      real(kind = 8) :: lower_encounter, upper_encounter
      real(kind = 8), dimension(125) :: spot_encounter, spot_predation
      real(kind = 8) :: encounter, r8_uniform_ab
      
      ! for linking the flounder to the spot
      integer :: i4_uniform_ab, lower_flounder_num, upper_flounder_num

      ! defining the variable for the daily simulations
      days = 20
      alive = 1
      lower_daily_mort = 0.000001
      upper_daily_mort = 1.000000
     
      ! defining variables for the truncated normal distribution for sizes
      n_spot = 125
      mu_small_spot = 42.8
      lower_small_spot = 33.0
      upper_small_spot = 49.0
      
      mu_large_spot = 63.7
      lower_large_spot = 54.0
      upper_large_spot = 78.0

      mu_small_flounder = 113.3
      lower_small_flounder = 92.0
      upper_small_flounder = 135.0

      mu_large_flounder = 181.5
      lower_large_flounder = 160.0
      upper_large_flounder = 200.0

      sigma_spot = 4.0
      sigma_flounder = 35.0
      seed = 123456789

      ! I have organized the different treatments as cases. The user will start by selecting the simulation that you want to run
      PRINT *, "Which model simulation would you like to run?" 
      PRINT *, "Type 1, 2, 3, or 4"
      PRINT *, "1 - Control"
      PRINT *, "2 - Small Flounder"
      PRINT *, "3 - Large Flounder" 
      PRINT *, "4 - Mixed"
      READ *, simulation
      
      ! Get the corresponding treatment
      SELECT CASE(simulation)
      
!----------------------------------------------------------------------------------------------------
!-------------------------------- Control Treatment: No flounder ------------------------------------
!----------------------------------------------------------------------------------------------------
         CASE(1)
         treatment = "Control, no Flounder"

         ! Size distributions
         ! SMALL SPOT
         write ( *, '(a)' ) ' '
         write ( *, '(a)' ) 'SIZES FOR SMALL SPOT '
         write ( *, '(a,g14.6)' ) '  Lower limit A =', lower_small_spot
         write ( *, '(a,g14.6)' ) '  Upper limit B =', upper_small_spot
         write ( *, '(a,g14.6)' ) '  MU =    ', mu_small_spot
         write ( *, '(a,g14.6)' ) '  SIGMA = ', sigma_spot
         write ( *, '(a,i12)' ) '  SEED = ', seed
         write ( *, '(a)' ) ''
         
         do i = 1,n_spot
           call truncated_normal(mu_small_spot,sigma_spot,
     &          lower_small_spot,upper_small_spot,seed,size_of_spot)
           WRITE (*,*) i, size_of_spot
           size_small_spot(i) = size_of_spot
         end do
           min = minval (size_small_spot, dim = 1)
           max = maxval (size_small_spot, dim = 1)
           WRITE(*,*) "Min Size = ", min
           WRITE(*,*) "Max Size = ", max
           
         ! LARGE SPOT 
         write ( *, '(a)' ) ' '
         write ( *, '(a,g14.6)' ) '  Lower limit A =', lower_large_spot
         write ( *, '(a,g14.6)' ) '  Upper limit B =', upper_large_spot
         write ( *, '(a,g14.6)' ) '  MU =    ', mu_large_spot
         write ( *, '(a,g14.6)' ) '  SIGMA = ', sigma_spot
         write ( *, '(a,i12)' ) '  SEED = ', seed
         write ( *, '(a)' ) ''
         
         do i = 1,n_spot
          call truncated_normal(mu_large_spot,sigma_spot,
     &         lower_large_spot,upper_large_spot,seed,size_of_spot)
          WRITE (*,*) i, size_of_spot
          size_large_spot(i) = size_of_spot
         end do
          min = minval (size_large_spot, dim = 1)
          max = maxval (size_large_spot, dim = 1)
          WRITE(*,*) ' '
          WRITE(*,*) "Min Size = ", min
          WRITE(*,*) "Max Size = ", max
           
         ! Expressing daily natural mortality of 0.018 for spot over the 20 day simulation. Have to randomly assign a worth value
         do i = 1, n_spot
            small_spot(4,i) = alive
            large_spot(4,i) = alive
         end do

         do i = 1, days
           do j = 1, n_spot
              small_spot(1,j) = j
              large_spot(1,j) = j
              small_spot(2,j) = size_small_spot(j)
              large_spot(2,j) = size_large_spot(j)
              small_spot(3,j) = r8_uniform_ab(lower_daily_mort, 
     &             upper_daily_mort, seed)
              large_spot(3,j) = r8_uniform_ab(lower_daily_mort, 
     &             upper_daily_mort, seed)
           if (small_spot(3,j) <= daily_mort) then
               small_spot(4,j) = 0
           end if
           if (large_spot(3,j) <= daily_mort) then
               large_spot(4,j) = 0
           end if
            end do
         end do
         
         WRITE(*,*) ' '
         WRITE(*,*) "The min values are", minval(small_spot, dim = 2)
         WRITE(*,*) "The max values are", maxval(small_spot, dim = 2)
         WRITE(*,*) ' '
         WRITE(*,*) "The min values are", minval(large_spot, dim = 2)
         WRITE(*,*) "The max values are", maxval(large_spot, dim = 2)

               
         ! Writing the final array size and checking to make sure that it matches the final value for n_spot
         WRITE(*,*) "The final number of spot is", n_spot  
         
         ! Writing the excel file
         OPEN(unit = 1, access = "sequential", action = "write", 
     &        status = "replace", 
     &        file = "Rice_model_control_small_spot.csv")
              WRITE(1,'(4(g14.7,:,","))') small_spot(1:4,:)
         OPEN(unit = 1, access = "sequential", action = "write", 
     &        status = "replace", 
     &        file = "Rice_model_control_large_spot.csv")
              WRITE(1,'(4(g14.7,:,","))') large_spot(1:4,:)
         CLOSE(1)

         WRITE(*,*) "The final number of large spot in case 1 is", 
     &        count(small_spot(4,:) .eq. 1)
         WRITE(*,*) "The final number of large spot in case 1 is", 
     &        count(large_spot(4,:) .eq. 1)

        

!----------------------------------------------------------------------------------------------------
!--------------------------- Small Flounder Treatment: 10 small flounder ----------------------------
!----------------------------------------------------------------------------------------------------
        CASE(2)
         treatment = "10 Small flounder"
         ! assigning values to the inputs, these numbers will change depending on the treatment
         num_of_small_flounder = 10
         num_of_large_flounder = 0
         lower_flounder_num = 1
         upper_flounder_num = 10
         
         ! the truncated normal distribution for randomly assigning sizes for small spot
         ! SMALL SPOT
         write ( *, '(a)' ) ' '
         write ( *, '(a)' ) 'SIZES FOR SMALL SPOT '
         write ( *, '(a,g14.6)' ) '  Lower limit A =', lower_small_spot
         write ( *, '(a,g14.6)' ) '  Upper limit B =', upper_small_spot
         write ( *, '(a,g14.6)' ) '  MU =    ', mu_small_spot
         write ( *, '(a,g14.6)' ) '  SIGMA = ', sigma_spot
         write ( *, '(a,i12)' ) '  SEED = ', seed
         write ( *, '(a)' ) ''        
         
         do i = 1,n_spot
           call truncated_normal(mu_small_spot,sigma_spot,
     &          lower_small_spot,upper_small_spot,seed,size_of_spot)
           WRITE (*,*) i, size_of_spot
           size_small_spot(i) = size_of_spot
         end do
           min = minval (size_small_spot, dim = 1)
           max = maxval (size_small_spot, dim = 1)
           WRITE(*,*) ' '
           WRITE(*,*) "Min Size = ", min
           WRITE(*,*) "Max Size = ", max

         ! LARGE SPOT
         write ( *, '(a)' ) ' '
         write ( *, '(a)' ) 'SIZES FOR LARGE SPOT '
         write ( *, '(a,g14.6)' ) '  Lower limit A =', lower_large_spot
         write ( *, '(a,g14.6)' ) '  Upper limit B =', upper_large_spot
         write ( *, '(a,g14.6)' ) '  MU =    ', mu_large_spot
         write ( *, '(a,g14.6)' ) '  SIGMA = ', sigma_spot
         write ( *, '(a,i12)' ) '  SEED = ', seed
         write ( *, '(a)' ) ''
      
         do i = 1,n_spot
          call truncated_normal(mu_large_spot,sigma_spot,
     &         lower_large_spot,upper_large_spot,seed,size_of_spot)
          WRITE (*,*) i, size_of_spot
          size_large_spot(i) = size_of_spot
         end do
          min = minval (size_large_spot, dim = 1)
          max = maxval (size_large_spot, dim = 1)
          WRITE (*,*) ' '
          WRITE(*,*) "Min Size = ", min
          WRITE(*,*) "Max Size = ", max

         ! SMALL FlOUNDER
         write ( *, '(a)' ) ' '
         write ( *, '(a)' ) 'SIZES FOR SMALL FLOUNDER '
         write ( *, '(a,g14.6)' )'Lower limit A =', lower_small_flounder
         write ( *, '(a,g14.6)' )'Upper limit B =', upper_small_flounder
         write ( *, '(a,g14.6)' )'MU =    ', mu_small_flounder
         write ( *, '(a,g14.6)' )'SIGMA = ', sigma_flounder
         write ( *, '(a,i12)' )'SEED = ', seed
         write ( *, '(a)' ) ''
         do i = 1,num_of_small_flounder
           call truncated_normal(mu_small_flounder,sigma_flounder,
     &          lower_small_flounder,upper_small_flounder,
     &          seed,size_of_flounder)
           WRITE (*,*) i, size_of_flounder
           small_flounder_case_2(i) = size_of_flounder
         end do
           min = minval (small_flounder_case_2, dim = 1)
           max = maxval (small_flounder_case_2, dim = 1)
           WRITE(*,*) ' '
           WRITE(*,*) "Min Size = ", min
           WRITE(*,*) "Max Size = ", max

         ! defining the variables for encounters
         lower_encounter = 0.000001
         upper_encounter = 1.000000

         ! creating the alive variable
         alive = 1
         do i = 1, n_spot
            small_spot(9,i) = alive
            large_spot(9,i) = alive
         end do

         ! Expressing daily natural mortality of 0.018 for spot over the 20 day simulation
         do i = 1, days
           do k = 1, num_of_small_flounder
              small_flounder(1,k) = k

              if (i == 1) then
              small_flounder(2,k) = small_flounder_case_2(k)
              else 
              small_flounder(2,k) = small_flounder(2,k) + 
     &                   growth_small_flounder
              end if
           end do
   
           do j = 1, n_spot
              small_spot(1,j) = j
              small_spot(2,j) = size_small_spot(j)
              large_spot(1,j) = j
              large_spot(2,j) = size_large_spot(j)

           ! if dead then nothing happens after this
           if (small_spot(9,j) == 1) then
           
           ! ENCOUNTER: assigning a number between 0 and 1 for the probability of an encounter, if the number is <= 0.07 for small spot, then an encounter occurred
              small_spot(3,j) = r8_uniform_ab(lower_encounter, 
     &             upper_encounter, seed)
              
              large_spot(3,j) = r8_uniform_ab(lower_encounter, 
     &             upper_encounter, seed)

           ! FLOUNDER NUMBER: assigning a rancom number between 0 and 9 to relate back to flounder
              small_spot(4,j) = i4_uniform_ab(lower_flounder_num, 
     &             upper_flounder_num)
              large_spot(4,j) = i4_uniform_ab(lower_flounder_num, 
     &             upper_flounder_num)
             
           ! MERGE WITH FLOUNDER
              if (small_spot(4,j) == 1) then
                 small_spot(5,j) = small_flounder(2,1)
              else if (small_spot(4,j) == 2) then
                 small_spot(5,j) = small_flounder(2,2)
              else if (small_spot(4,j) == 3) then
                 small_spot(5,j) = small_flounder(2,3)
              else if (small_spot(4,j) == 4) then
                 small_spot(5,j) = small_flounder(2,4)
              else if (small_spot(4,j) == 5) then
                 small_spot(5,j) = small_flounder(2,5)
              else if (small_spot(4,j) == 6) then
                 small_spot(5,j) = small_flounder(2,6)
              else if (small_spot(4,j) == 7) then
                 small_spot(5,j) = small_flounder(2,7)
              else if (small_spot(4,j) == 8) then
                 small_spot(5,j) = small_flounder(2,8)
              else if (small_spot(4,j) == 9) then
                 small_spot(5,j) = small_flounder(2,9)
              else if (small_spot(4,j) == 10) then
                 small_spot(5,j) = small_flounder(2,10)
              end if

              if (large_spot(4,j) == 1) then
                 large_spot(5,j) = small_flounder(2,1)
              else if (large_spot(4,j) == 2) then
                 large_spot(5,j) = small_flounder(2,2)
              else if (large_spot(4,j) == 3) then
                 large_spot(5,j) = small_flounder(2,3)
              else if (large_spot(4,j) == 4) then
                 large_spot(5,j) = small_flounder(2,4)
              else if (large_spot(4,j) == 5) then
                 large_spot(5,j) = small_flounder(2,5)
              else if (large_spot(4,j) == 6) then
                 large_spot(5,j) = small_flounder(2,6)
              else if (large_spot(4,j) == 7) then
                 large_spot(5,j) = small_flounder(2,7)
              else if (large_spot(4,j) == 8) then
                 large_spot(5,j) = small_flounder(2,8)
              else if (large_spot(4,j) == 9) then
                 large_spot(5,j) = small_flounder(2,9)
              else if (large_spot(4,j) == 10) then
                 large_spot(5,j) = small_flounder(2,10)
              end if

           ! PROFITABILITY RATIO
              small_spot(6,j) = small_spot(5,j) /
     &             small_spot(2,j)
              large_spot(6,j) = large_spot(5,j) /
     &             large_spot(2,j)

           ! PREDATION
              small_spot(7,j) = r8_uniform_ab(lower_encounter,
     &             upper_encounter, seed)
              large_spot(7,j) = r8_uniform_ab(lower_encounter, 
     &             upper_encounter, seed)

           ! NATURAL MORTALITY - this shouldn't be random. it's a set value
              small_spot(8,j) = r8_uniform_ab(lower_daily_mort, 
     &             upper_daily_mort, seed)
              large_spot(8,j) = r8_uniform_ab(lower_daily_mort, 
     &             upper_daily_mort, seed)
          
           ! ALIVE = 0 STATEMENTS
              if (small_spot(8,j) <= daily_mort .or. small_spot(3,j) 
     &          <= prob_of_enc_small .and. small_spot(6,j) >= ratio_min 
     &          .and. small_spot(7,j) > prob_predation) then
                small_spot(9,j) = 0
              end if

              if (large_spot(8,j) <= daily_mort .or. large_spot(3,j) 
     &             <= prob_of_enc_large.and.large_spot(6,j) >= ratio_min 
     &             .and. large_spot(7,j) > prob_predation) then
                large_spot(9,j) = 0
              end if
             end if
           end do 
            WRITE(*,*) ' '
            WRITE(*,*) "The min values are", minval(small_spot, dim = 2)
            WRITE(*,*) "The max values are", maxval(small_spot, dim = 2)
            WRITE(*,*) ' '
            WRITE(*,*) "The min values are", minval(large_spot, dim = 2)
            WRITE(*,*) "The max values are", maxval(large_spot, dim = 2)
      
         ! Writing the excel file for small and large spot
         OPEN(unit = 1, access = "sequential", action = "write", 
     &        status = "replace", 
     &        file = "Rice_model_case_2_small_spot.csv")
              WRITE(1,'(9(g14.7,:,","))') small_spot
         CLOSE(1)         
         OPEN(unit = 1, access = "sequential", action = "write", 
     &        status = "replace", 
     &        file = "Rice_model_case_2_large_spot.csv")
              WRITE(1,'(9(g14.7,:,","))') large_spot
         CLOSE(1)
         end do
         
         !PRINT *, ((small_spot(i,j),j=1,n_spot), i = 1,3OA)
         !PRINT *, ((large_spot(i,j),j=1,n_spot), i = 1,3OA)         
         

               
         ! Writing the final array size and checking to make sure that it matches the final value for n_spot
         WRITE(*,*) "The final number of small spot in case 2 is", 
     &        count(small_spot(9,:) .eq. 1)
         WRITE(*,*) "The final number of large spot in case 2 is", 
     &        count(large_spot(9,:) .eq. 1)

        


!----------------------------------------------------------------------------------------------------
!--------------------------- Large Flounder Treatment: 4 large flounder -----------------------------
!----------------------------------------------------------------------------------------------------
        CASE(3)
         treatment = "4 Large flounder"
         ! assigning values to the inputs, these numbers will change depending on the treatment
         num_of_small_flounder = 0
         num_of_large_flounder = 4
         lower_flounder_num = 1
         upper_flounder_num = 4
         ! the truncated normal ditribution for randomaly assigning sizes for small spot
         write ( *, '(a)' ) ' '
         write ( *, '(a,g14.6)' ) '  Lower limit A =    ', a
         write ( *, '(a,g14.6)' ) '  Upper limit B =    ', b
         write ( *, '(a,g14.6)' ) '  MU =    ', mu
         write ( *, '(a,g14.6)' ) '  SIGMA = ', sigma
         write ( *, '(a,i12)' ) '  SEED = ', seed
         write ( *, '(a)' ) ''
   
         
         ! SMALL SPOT
         do i = 1,n_spot
           call truncated_normal(mu_small_spot,sigma_spot,
     &         lower_small_spot,upper_small_spot,seed,size_of_spot)
           WRITE (*,*) i, size_of_spot
           size_small_spot(i) = size_of_spot
         end do
           min = minval (size_small_spot, dim = 1)
           max = maxval (size_small_spot, dim = 1)
           WRITE(*,*) ' '
           WRITE(*,*) "Min Size = ", min
           WRITE(*,*) "Max Size = ", max

         ! LARGE SPOT
         do i = 1,n_spot
          call truncated_normal(mu_large_spot,sigma_spot,
     &           lower_large_spot,upper_large_spot,seed,size_of_spot)
          WRITE (*,*) i, size_of_spot
          size_large_spot(i) = size_of_spot
         end do
          min = minval (size_large_spot, dim = 1)
          max = maxval (size_large_spot, dim = 1)
          WRITE(*,*) ' '
          WRITE(*,*) "Min Size = ", min
          WRITE(*,*) "Max Size = ", max

         ! LARGE FlOUNDER
         write ( *, '(a)' ) ' '
         write ( *, '(a)' ) 'SIZES FOR LARGE FLOUNDER '
         write ( *, '(a,g14.6)' )'Lower limit A =', lower_large_flounder
         write ( *, '(a,g14.6)' )'Upper limit B =', upper_large_flounder
         write ( *, '(a,g14.6)' )'MU =    ', mu_large_flounder
         write ( *, '(a,g14.6)' )'SIGMA = ', sigma_flounder
         write ( *, '(a,i12)' )'SEED = ', seed
         write ( *, '(a)' ) ''
         do i = 1,num_of_large_flounder
           call truncated_normal(mu_large_flounder,sigma_flounder,
     &         lower_large_flounder,upper_large_flounder,
     &         seed,size_of_flounder)
           WRITE (*,*) i, size_of_flounder
           large_flounder_case_3(i) = size_of_flounder
         end do
           min = minval (large_flounder_case_3, dim = 1)
           max = maxval (large_flounder_case_3, dim = 1)
           WRITE(*,*) ' '
           WRITE(*,*) "Min Size = ", min
           WRITE(*,*) "Max Size = ", max

         ! defining the variables for encounters
         lower_encounter = 0.000001
         upper_encounter = 1.000000

         ! creating the alive variable
         alive = 1
         do i = 1, n_spot
            small_spot(9,i) = alive
            large_spot(9,i) = alive
         end do

         ! 20 day simulation
         do i = 1, days
           do k = 1, num_of_large_flounder
              large_flounder(1,k) = k

              if (i == 1) then
              large_flounder(2,k) = large_flounder_case_3(k)
              else 
              large_flounder(2,k) = large_flounder(2,k) + 
     &                   growth_large_flounder
              end if
           end do
   
           do j = 1, n_spot
              small_spot(1,j) = j
              small_spot(2,j) = size_small_spot(j)
              large_spot(1,j) = j
              large_spot(2,j) = size_large_spot(j)

           ! if dead then nothing happens after this
           if (small_spot(9,j) == 1) then
           
           ! ENCOUNTER: assigning a number between 0 and 1 for the probability of an encounter, if the number is <= 0.07 for small spot, then an encounter occurred
              small_spot(3,j) = r8_uniform_ab(lower_encounter, 
     &             upper_encounter, seed)
              
              large_spot(3,j) = r8_uniform_ab(lower_encounter, 
     &             upper_encounter, seed)

           ! FLOUNDER NUMBER: assigning a rancom number between 0 and 9 to relate back to flounder
              small_spot(4,j) = i4_uniform_ab(lower_flounder_num, 
     &             upper_flounder_num)
              large_spot(4,j) = i4_uniform_ab(lower_flounder_num, 
     &             upper_flounder_num)
             
           ! MERGE WITH FLOUNDER
              if (small_spot(4,j) == 1) then
                 small_spot(5,j) = large_flounder(2,1)
              else if (small_spot(4,j) == 2) then
                 small_spot(5,j) = large_flounder(2,2)
              else if (small_spot(4,j) == 3) then
                 small_spot(5,j) = large_flounder(2,3)
              else if (small_spot(4,j) == 4) then
                 small_spot(5,j) = large_flounder(2,4)
              end if

              if (large_spot(4,j) == 1) then
                 large_spot(5,j) = large_flounder(2,1)
              else if (large_spot(4,j) == 2) then
                 large_spot(5,j) = large_flounder(2,2)
              else if (large_spot(4,j) == 3) then
                 large_spot(5,j) = large_flounder(2,3)
              else if (large_spot(4,j) == 4) then
                 large_spot(5,j) = large_flounder(2,4)
              end if

           ! PROFITABILITY RATIO
              small_spot(6,j) = small_spot(5,j) /
     &             small_spot(2,j)
              large_spot(6,j) = large_spot(5,j) /
     &             large_spot(2,j)

           ! PREDATION
              small_spot(7,j) = r8_uniform_ab(lower_encounter,
     &             upper_encounter, seed)
              large_spot(7,j) = r8_uniform_ab(lower_encounter, 
     &             upper_encounter, seed)

           ! NATURAL MORTALITY - this shouldn't be random. it's a set value
              small_spot(8,j) = r8_uniform_ab(lower_daily_mort, 
     &             upper_daily_mort, seed)
              large_spot(8,j) = r8_uniform_ab(lower_daily_mort, 
     &             upper_daily_mort, seed)
          
           ! ALIVE = 0 STATEMENTS
              if (small_spot(8,j) <= daily_mort .or. small_spot(3,j) 
     &          <= prob_of_enc_small .and. small_spot(6,j) >= ratio_min 
     &          .and. small_spot(7,j) > prob_predation) then
                small_spot(9,j) = 0
              end if

              if (large_spot(8,j) <= daily_mort .or. large_spot(3,j) 
     &             <= prob_of_enc_large.and.large_spot(6,j) >= ratio_min 
     &             .and. large_spot(7,j) > prob_predation) then
                large_spot(9,j) = 0
              end if
             end if
           end do 

            WRITE(*,*) ' '
            WRITE(*,*) "The min values are", minval(small_spot, dim = 2)
            WRITE(*,*) "The max values are", maxval(small_spot, dim = 2)
            WRITE(*,*) ' '
            WRITE(*,*) "The min values are", minval(large_spot, dim = 2)
            WRITE(*,*) "The max values are", maxval(large_spot, dim = 2)
      
         ! Writing the excel file for small and large spot
         OPEN(unit = 1, access = "sequential", action = "write", 
     &        status = "replace", 
     &        file = "Rice_model_case_3_small_spot.csv")
              WRITE(1,'(9(g14.7,:,","))') small_spot
         CLOSE(1)         
         OPEN(unit = 1, access = "sequential", action = "write", 
     &        status = "replace", 
     &        file = "Rice_model_case_3_large_spot.csv")
              WRITE(1,'(9(g14.7,:,","))') large_spot
         CLOSE(1)
         end do
         
         !PRINT *, ((small_spot(i,j),j=1,n_spot), i = 1,3OA)
         !PRINT *, ((large_spot(i,j),j=1,n_spot), i = 1,3OA)         
         

               
         ! Writing the final array size and checking to make sure that it matches the final value for n_spot
         WRITE(*,*) "The final number of small spot in case 3 is", 
     &        count(small_spot(9,:) .eq. 1)
         WRITE(*,*) "The final number of large spot in case 3 is", 
     &        count(large_spot(9,:) .eq. 1)

         

!----------------------------------------------------------------------------------------------------
!--------------------------- Mixed Treatment: 5 small and 2 large flounder --------------------------
!----------------------------------------------------------------------------------------------------
        CASE(4) 
         treatment = "Mixed, 5 small and 2 large flounder"
         ! assigning values to the inputs, these numbers will change depending on the treatment
         num_of_small_flounder = 5
         num_of_large_flounder = 2
         lower_flounder_num = 1
         upper_flounder_num = 7
         ! the truncated normal ditribution for randomaly assigning sizes for small spot
         write ( *, '(a)' ) ' '
         write ( *, '(a,g14.6)' ) '  Lower limit A =    ', a
         write ( *, '(a,g14.6)' ) '  Upper limit B =    ', b
         write ( *, '(a,g14.6)' ) '  MU =    ', mu
         write ( *, '(a,g14.6)' ) '  SIGMA = ', sigma
         write ( *, '(a,i12)' ) '  SEED = ', seed
         write ( *, '(a)' ) ''
   
         
         ! SMALL SPOT
         do i = 1,n_spot
           call truncated_normal(mu_small_spot,sigma_spot,
     &         lower_small_spot,upper_small_spot,seed,size_of_spot)
           WRITE (*,*) i, size_of_spot
           size_small_spot(i) = size_of_spot
         end do
           min = minval (size_small_spot, dim = 1)
           max = maxval (size_small_spot, dim = 1)
           WRITE(*,*) ' '
           WRITE(*,*) "Min Size = ", min
           WRITE(*,*) "Max Size = ", max

         ! LARGE SPOT
         do i = 1,n_spot
          call truncated_normal(mu_large_spot,sigma_spot,
     &           lower_large_spot,upper_large_spot,seed,size_of_spot)
          WRITE (*,*) i, size_of_spot
          size_large_spot(i) = size_of_spot
         end do
          min = minval (size_large_spot, dim = 1)
          max = maxval (size_large_spot, dim = 1)
          WRITE(*,*) ' '
          WRITE(*,*) "Min Size = ", min
          WRITE(*,*) "Max Size = ", max

         ! SMALL FlOUNDER
         write ( *, '(a)' ) ' '
         write ( *, '(a)' ) 'SIZES FOR SMALL FLOUNDER '
         write ( *, '(a,g14.6)' )'Lower limit A =', lower_small_flounder
         write ( *, '(a,g14.6)' )'Upper limit B =', upper_small_flounder
         write ( *, '(a,g14.6)' )'MU =    ', mu_small_flounder
         write ( *, '(a,g14.6)' )'SIGMA = ', sigma_flounder
         write ( *, '(a,i12)' )'SEED = ', seed
         write ( *, '(a)' ) ''
         do i = 1,num_of_small_flounder
           call truncated_normal(mu_small_flounder,sigma_flounder,
     &         lower_small_flounder,upper_small_flounder,
     &         seed,size_of_flounder)
           WRITE (*,*) i, size_of_flounder
           small_flounder_case_4(i) = size_of_flounder
         end do
           min = minval (small_flounder_case_4, dim = 1)
           max = maxval (small_flounder_case_4, dim = 1)
           WRITE(*,*) ' '
           WRITE(*,*) "Min Size = ", min
           WRITE(*,*) "Max Size = ", max

         ! LARGE FlOUNDER
         write ( *, '(a)' ) ' '
         write ( *, '(a)' ) 'SIZES FOR LARGE FLOUNDER '
         write ( *, '(a,g14.6)' )'Lower limit A =', lower_large_flounder
         write ( *, '(a,g14.6)' )'Upper limit B =', upper_large_flounder
         write ( *, '(a,g14.6)' )'MU =    ', mu_large_flounder
         write ( *, '(a,g14.6)' )'SIGMA = ', sigma_flounder
         write ( *, '(a,i12)' )'SEED = ', seed
         write ( *, '(a)' ) ''
         do i = 1,num_of_large_flounder
           call truncated_normal(mu_large_flounder,sigma_flounder,
     &         lower_large_flounder,upper_large_flounder,
     &         seed,size_of_flounder)
           WRITE (*,*) i, size_of_flounder
           large_flounder_case_4(i) = size_of_flounder
         end do
           min = minval (large_flounder_case_4, dim = 1)
           max = maxval (large_flounder_case_4, dim = 1)
           WRITE(*,*) ' '
           WRITE(*,*) "Min Size = ", min
           WRITE(*,*) "Max Size = ", max

         ! defining the variables for encounters
         lower_encounter = 0.000001
         upper_encounter = 1.000000

         ! creating the alive variable
         alive = 1
         do i = 1, n_spot
            small_spot(9,i) = alive
            large_spot(9,i) = alive
         end do

         ! Expressing daily natural mortality of 0.018 for spot over the 20 day simulation
         do i = 1, days
           ! DAILY GROWTH FOR SMALL AND LARGE FLOUNDER
           do k = 1, num_of_small_flounder
              small_flounder(1,k) = k

              if (i == 1) then
              small_flounder(2,k) = small_flounder_case_4(k)
               else 
              small_flounder(2,k) = small_flounder(2,k) + 
     &                   growth_small_flounder
              end if
           end do
           do k = 1, num_of_large_flounder
              large_flounder(1,k) = k

              if (i == 1) then
              large_flounder(2,k) = large_flounder_case_4(k)
               else 
              large_flounder(2,k) = large_flounder(2,k) + 
     &                   growth_large_flounder
              end if
           end do   
           do j = 1, n_spot
              small_spot(1,j) = j
              small_spot(2,j) = size_small_spot(j)
              large_spot(1,j) = j
              large_spot(2,j) = size_large_spot(j)

           ! if dead then nothing happens after this
           if (small_spot(9,j) == 1) then
           
           ! ENCOUNTER: assigning a number between 0 and 1 for the probability of an encounter, if the number is <= 0.07 for small spot, then an encounter occurred
              small_spot(3,j) = r8_uniform_ab(lower_encounter, 
     &             upper_encounter, seed)
              
              large_spot(3,j) = r8_uniform_ab(lower_encounter, 
     &             upper_encounter, seed)

           ! FLOUNDER NUMBER: assigning a rancom number between 0 and 9 to relate back to flounder
              small_spot(4,j) = i4_uniform_ab(lower_flounder_num, 
     &             upper_flounder_num)
              large_spot(4,j) = i4_uniform_ab(lower_flounder_num, 
     &             upper_flounder_num)
             
           ! MERGE WITH FLOUNDER
              if (small_spot(4,j) == 1) then
                 small_spot(5,j) = small_flounder(2,1)
              else if (small_spot(4,j) == 2) then
                 small_spot(5,j) = small_flounder(2,2)
              else if (small_spot(4,j) == 3) then
                 small_spot(5,j) = small_flounder(2,3)
              else if (small_spot(4,j) == 4) then
                 small_spot(5,j) = small_flounder(2,4)
              else if (small_spot(4,j) == 5) then
                 small_spot(5,j) = small_flounder(2,5)
              else if (small_spot(4,j) == 6) then
                 small_spot(5,j) = large_flounder(2,1)
              else if (small_spot(4,j) == 7) then
                 small_spot(5,j) = large_flounder(2,2)
              end if

              if (large_spot(4,j) == 1) then
                 large_spot(5,j) = small_flounder(2,1)
              else if (large_spot(4,j) == 2) then
                 large_spot(5,j) = small_flounder(2,2)
              else if (large_spot(4,j) == 3) then
                 large_spot(5,j) = small_flounder(2,3)
              else if (large_spot(4,j) == 4) then
                 large_spot(5,j) = small_flounder(2,4)
              else if (large_spot(4,j) == 5) then
                 large_spot(5,j) = small_flounder(2,5)
              else if (large_spot(4,j) == 6) then
                 large_spot(5,j) = large_flounder(2,1)
              else if (large_spot(4,j) == 7) then
                 large_spot(5,j) = large_flounder(2,2)
              end if

           ! PROFITABILITY RATIO
              small_spot(6,j) = small_spot(5,j) /
     &             small_spot(2,j)
              large_spot(6,j) = large_spot(5,j) /
     &             large_spot(2,j)

           ! PREDATION
              small_spot(7,j) = r8_uniform_ab(lower_encounter,
     &             upper_encounter, seed)
              large_spot(7,j) = r8_uniform_ab(lower_encounter, 
     &             upper_encounter, seed)

           ! NATURAL MORTALITY - this shouldn't be random. it's a set value
              small_spot(8,j) = r8_uniform_ab(lower_daily_mort, 
     &             upper_daily_mort, seed)
              large_spot(8,j) = r8_uniform_ab(lower_daily_mort, 
     &             upper_daily_mort, seed)
          
           ! ALIVE = 0 STATEMENTS
              if (small_spot(8,j) <= daily_mort .or. small_spot(3,j) 
     &          <= prob_of_enc_small .and. small_spot(6,j) >= ratio_min 
     &          .and. small_spot(7,j) > prob_predation) then
                small_spot(9,j) = 0
              end if

              if (large_spot(8,j) <= daily_mort .or. large_spot(3,j) 
     &             <= prob_of_enc_large.and.large_spot(6,j) >= ratio_min 
     &             .and. large_spot(7,j) > prob_predation) then
                large_spot(9,j) = 0
              end if
             end if
           end do 
            WRITE(*,*) ' '
            WRITE(*,*) "The min values are", minval(small_spot, dim = 2)
            WRITE(*,*) "The max values are", maxval(small_spot, dim = 2)
            WRITE(*,*) ' '
            WRITE(*,*) "The min values are", minval(large_spot, dim = 2)
            WRITE(*,*) "The max values are", maxval(large_spot, dim = 2)
      
         ! Writing the excel file for small and large spot
         OPEN(unit = 1, access = "sequential", action = "write", 
     &        status = "replace", 
     &        file = "Rice_model_case_4_small_spot.csv")
              WRITE(1,'(9(g14.7,:,","))') small_spot
         CLOSE(1)         
         OPEN(unit = 1, access = "sequential", action = "write", 
     &        status = "replace", 
     &        file = "Rice_model_case_4_large_spot.csv")
              WRITE(1,'(9(g14.7,:,","))') large_spot
         CLOSE(1)
         end do
         
         !PRINT *, ((small_spot(i,j),j=1,n_spot), i = 1,3OA)
         !PRINT *, ((large_spot(i,j),j=1,n_spot), i = 1,3OA)         
         

               
         ! Writing the final array size and checking to make sure that it matches the final value for n_spot
         WRITE(*,*) "The final number of small spot in case 4 is", 
     &        count(small_spot(9,:) .eq. 1)
         WRITE(*,*) "The final number of large spot in case 4 is", 
     &        count(large_spot(9,:) .eq. 1)
      END SELECT
      WRITE (*,*) "Treatment = ",treatment
         

      END PROGRAM Rice

         !creating the truncated normal distribution for both the small and large spot
         SUBROUTINE truncated_normal (mu,sigma,a,b,seed,x)
          IMPLICIT NONE
          real(kind = 8) a, b, mu, sigma, x
          real(kind = 8) alpha,alpha_cdf,beta,beta_cdf,r8_uniform_01
          real(kind = 8) u,size_of_spot,xi,xi_cdf
          integer(kind = 4) :: seed
          
          alpha = (a - mu) / sigma
          beta = (b - mu) / sigma
        
          call normal_01_cdf (alpha,alpha_cdf)
          call normal_01_cdf (beta,beta_cdf)
          
          u = r8_uniform_01 (seed)
          xi_cdf = alpha_cdf + u * (beta_cdf - alpha_cdf)
          call normal_01_cdf_inv (xi_cdf, xi)

          x = mu + sigma * xi
   
         return
         END SUBROUTINE truncated_normal
      
        subroutine normal_01_cdf (x, cdf)
         implicit none

         real ( kind = 8 ), parameter :: a1 = 0.398942280444D+00
         real ( kind = 8 ), parameter :: a2 = 0.399903438504D+00
         real ( kind = 8 ), parameter :: a3 = 5.75885480458D+00
         real ( kind = 8 ), parameter :: a4 = 29.8213557808D+00
         real ( kind = 8 ), parameter :: a5 = 2.62433121679D+00
         real ( kind = 8 ), parameter :: a6 = 48.6959930692D+00
         real ( kind = 8 ), parameter :: a7 = 5.92885724438D+00
         real ( kind = 8 ), parameter :: b0 = 0.398942280385D+00
         real ( kind = 8 ), parameter :: b1 = 3.8052D-08
         real ( kind = 8 ), parameter :: b2 = 1.00000615302D+00
         real ( kind = 8 ), parameter :: b3 = 3.98064794D-04
         real ( kind = 8 ), parameter :: b4 = 1.98615381364D+00
         real ( kind = 8 ), parameter :: b5 = 0.151679116635D+00
         real ( kind = 8 ), parameter :: b6 = 5.29330324926D+00
         real ( kind = 8 ), parameter :: b7 = 4.8385912808D+00
         real ( kind = 8 ), parameter :: b8 = 15.1508972451D+00
         real ( kind = 8 ), parameter :: b9 = 0.742380924027D+00
         real ( kind = 8 ), parameter :: b10 = 30.789933034D+00
         real ( kind = 8 ), parameter :: b11 = 3.99019417011D+00
         real ( kind = 8 ) cdf
         real ( kind = 8 ) q
         real ( kind = 8 ) x
         real ( kind = 8 ) y
         !  |X| <= 1.28.
    
      if ( abs ( x ) <= 1.28D+00 ) then
            
      y = 0.5D+00 * x * x
            
      q = 0.5D+00-abs(x)*(a1-a2*y/(y+a3-a4/(y+a5+a6/(y + a7))))
!
!  1.28 < |X| <= 12.7
!
      else if ( abs ( x ) <= 12.7D+00 ) then

      y = 0.5D+00 * x * x

      q=exp(-y)*b0/(abs(x)-b1+b2/(abs(x)+b3+b4/(abs(x)-b5+b6/(abs(x)+b7- 
     &  b8 / ( abs ( x ) + b9 + b10 / ( abs ( x ) + b11 ) ) ) ) ) )
!
!  12.7 < |X|
!
      else

      q = 0.0D+00

      end if
!
!  Take account of negative X.
!
      if ( x < 0.0D+00 ) then
        cdf = q
      else
        cdf = 1.0D+00 - q
      end if

      return
      END SUBROUTINE normal_01_cdf

      SUBROUTINE normal_01_cdf_inv (p,x)
      implicit none

      real ( kind = 8 ), parameter, dimension ( 8 ) :: a = (/ 
     &  3.3871328727963666080D+00, 
     &  1.3314166789178437745D+02, 
     &  1.9715909503065514427D+03, 
     &  1.3731693765509461125D+04, 
     &  4.5921953931549871457D+04, 
     &  6.7265770927008700853D+04, 
     &  3.3430575583588128105D+04, 
     &  2.5090809287301226727D+03 /)
      real ( kind = 8 ), parameter, dimension ( 8 ) :: b = (/ 
     &  1.0D+00, 
     &  4.2313330701600911252D+01, 
     &  6.8718700749205790830D+02, 
     &  5.3941960214247511077D+03, 
     &  2.1213794301586595867D+04, 
     &  3.9307895800092710610D+04, 
     &  2.8729085735721942674D+04, 
     &  5.2264952788528545610D+03 /)
      real ( kind = 8 ), parameter, dimension ( 8 ) :: c = (/ 
     &  1.42343711074968357734D+00, 
     &  4.63033784615654529590D+00, 
     &  5.76949722146069140550D+00, 
     &  3.64784832476320460504D+00, 
     &  1.27045825245236838258D+00, 
     &  2.41780725177450611770D-01, 
     &  2.27238449892691845833D-02, 
     &  7.74545014278341407640D-04 /)
      real ( kind = 8 ), parameter :: const1 = 0.180625D+00
      real ( kind = 8 ), parameter :: const2 = 1.6D+00
      real ( kind = 8 ), parameter, dimension ( 8 ) :: d = (/ 
     &  1.0D+00, 
     &  2.05319162663775882187D+00, 
     &  1.67638483018380384940D+00, 
     &  6.89767334985100004550D-01, 
     &  1.48103976427480074590D-01, 
     &  1.51986665636164571966D-02, 
     &  5.47593808499534494600D-04, 
     &  1.05075007164441684324D-09 /)
      real ( kind = 8 ), parameter, dimension ( 8 ) :: e = (/ 
     &  6.65790464350110377720D+00, 
     &  5.46378491116411436990D+00, 
     &  1.78482653991729133580D+00, 
     &  2.96560571828504891230D-01, 
     &  2.65321895265761230930D-02, 
     &  1.24266094738807843860D-03, 
     &  2.71155556874348757815D-05, 
     &  2.01033439929228813265D-07 /)
      real ( kind = 8 ), parameter, dimension ( 8 ) :: f = (/ 
     &  1.0D+00, 
     &  5.99832206555887937690D-01, 
     &  1.36929880922735805310D-01, 
     &  1.48753612908506148525D-02, 
     &  7.86869131145613259100D-04, 
     &  1.84631831751005468180D-05, 
     &  1.42151175831644588870D-07, 
     &  2.04426310338993978564D-15 /)
      real ( kind = 8 ) p
      real ( kind = 8 ) q
      real ( kind = 8 ) r
      real ( kind = 8 ) r8poly_value_horner
      real ( kind = 8 ), parameter :: split1 = 0.425D+00
      real ( kind = 8 ), parameter :: split2 = 5.0D+00
      real ( kind = 8 ) x

      if ( p <= 0.0D+00 ) then
      x = - huge ( x )
      return
      end if

      if ( 1.0D+00 <= p ) then
      x = huge ( x )
      return
      end if

      q = p - 0.5D+00

      if ( abs ( q ) <= split1 ) then

      r = const1 - q * q
      x = q * r8poly_value_horner ( 7, a, r ) 
     &     / r8poly_value_horner ( 7, b, r )

       else

      if ( q < 0.0D+00 ) then
      r = p
      else
      r = 1.0D+00 - p
      end if

      if ( r <= 0.0D+00 ) then

      x = huge ( x )

      else

      r = sqrt ( - log ( r ) )

      if ( r <= split2 ) then

        r = r - const2
        x = r8poly_value_horner ( 7, c, r ) 
     &     / r8poly_value_horner ( 7, d, r )

      else

        r = r - split2
        x = r8poly_value_horner ( 7, e, r ) 
     &     / r8poly_value_horner ( 7, f, r )

      end if

      end if

      if ( q < 0.0D+00 ) then
      x = -x
      end if

      end if

      return
      END SUBROUTINE normal_01_cdf_inv
      
      SUBROUTINE r8_uniform_ab_test
      implicit none
        
      integer, parameter :: rk = kind ( 1.0D+00 )
      
      real ( kind = rk ) a
      real ( kind = rk ) b
      real ( kind = rk ) r8_uniform_ab
      integer i
      integer seed
      
      a = 5.0D+00
      b = 10.0D+00
      seed = 123456789
      
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8_UNIFORM_AB_TEST'
      write ( *, '(a)' ) '  R8_UNIFORM_AB computes pseudorandom values '
      write ( *, '(a)' ) '  in an interval [A,B].'
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  The lower endpoint A = ', a
      write ( *, '(a,g14.6)' ) '  The upper endpoint B = ', b
      write ( *, '(a,i12)' ) '  The initial seed is ', seed
      
      write ( *, '(a)' ) ' '
      do i = 1, 10
      write ( *, '(2x,i8,2x,g14.6)' ) i, r8_uniform_ab ( a, b, seed )
      end do
      
      return
      end subroutine r8_uniform_ab_test


      function r8poly_value_horner(m,c,x)
        implicit none

        integer ( kind = 4 ) m
      
        real ( kind = 8 ) c(0:m)
        integer ( kind = 4 ) i
        real ( kind = 8 ) r8poly_value_horner
        real ( kind = 8 ) value
        real ( kind = 8 ) x

        value = c(m)
       do i = m - 1, 0, -1
        value = value * x + c(i)
       end do

       r8poly_value_horner = value

      return
      END FUNCTION r8poly_value_horner


      FUNCTION r8_uniform_01(seed)
      implicit none
        
      integer, parameter :: i4_huge = 2147483647
      integer k
      real (kind = 8) r8_uniform_01
      integer seed
      
      if ( seed == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8_UNIFORM_01 - Fatal error!'
      write ( *, '(a)' ) '  Input value of SEED = 0.'
      stop 1
      end if
      
      k = seed / 127773
      
      seed = 16807 * ( seed - k * 127773 ) - k * 2836
      
      if ( seed < 0 ) then
      seed = seed + i4_huge
      end if
      
      r8_uniform_01 = real (seed, kind = 8) * 4.656612875D-10
      
      return
      END FUNCTION r8_uniform_01

      function r8_uniform_ab(a,b,seed)
      implicit none

      real (kind = 8) a
      real (kind = 8) b
      integer, parameter :: i4_huge = 2147483647
      integer k
      real (kind = 8) r8_uniform_ab
      integer seed

      if ( seed == 0 ) then
         write ( *, '(a)' ) ' '
         write ( *, '(a)' ) 'R8_UNIFORM_AB - Fatal error!'
         write ( *, '(a)' ) '  Input value of SEED = 0.'
         stop 1
      end if

      k = seed / 127773
      
      seed = 16807 * ( seed - k * 127773 ) - k * 2836
      
      if ( seed < 0 ) then
         seed = seed + i4_huge
      end if
      
      r8_uniform_ab = a + ( b - a ) * 
     &   real (seed, kind = 8) * 4.656612875D-10
      
      return
      END FUNCTION r8_uniform_ab

      function i4_uniform_ab(a,b)
      implicit none
      integer  a
      integer  b
      integer  i4_uniform_ab
      real (kind = 8) r

      call random_number(harvest = r)
      i4_uniform_ab = a + int((b + 1 - a) * r)
      
      return
      END FUNCTION i4_uniform_ab

