           call product(ur,ui,that,uptr,upti,it)
c
c-----------------------------------------------------------------------------------------
c
c  create staple in the positive x direction
c
           tmpr = 0.0d0
           tmpi = 0.0d0

           ustr = cshift(uptr,dim=xdir,shift= ixq)
           usti = cshift(upti,dim=xdir,shift= ixq)

           do ic = 1,nc
              do jc = 1,nc
                 do kc = 1,nc

                    tmpr(:,:,:,:,ic,jc) = tmpr(:,:,:,:,ic,jc)
     &                   + upxr(:,:,:,:,ic,kc) * ustr(:,:,:,:,kc,jc) - upxi(:,:,:,:,ic,kc) * usti(:,:,:,:,kc,jc)

                    tmpi(:,:,:,:,ic,jc) = tmpi(:,:,:,:,ic,jc)
     &                   + upxr(:,:,:,:,ic,kc) * usti(:,:,:,:,kc,jc) + upxi(:,:,:,:,ic,kc) * ustr(:,:,:,:,kc,jc)

                 end do
              end do
           end do
c
c     go backwards from the positive x direction, forming staple
c
           sxpr = 0.0d0
           sxpi = 0.0d0

           usr  = cshift(upxr,dim=that,shift= it)
           usi  = cshift(upxi,dim=that,shift= it)

           do ic = 1,nc
              do jc = 1,nc
                 do kc = 1,nc

                    sxpr(:,:,:,:,ic,jc) = sxpr(:,:,:,:,ic,jc)
     &                   + tmpr(:,:,:,:,ic,kc) * usr(:,:,:,:,jc,kc) + tmpi(:,:,:,:,ic,kc) * usi(:,:,:,:,jc,kc)

                    sxpi(:,:,:,:,ic,jc) = sxpi(:,:,:,:,ic,jc)
     &                   - tmpr(:,:,:,:,ic,kc) * usi(:,:,:,:,jc,kc) + tmpi(:,:,:,:,ic,kc) * usr(:,:,:,:,jc,kc)

                 end do
              end do
           end do
c
c---------------------------------------------------------------------------
c
c  create staple in the negative x direction
c
           tmpr = 0.0d0
           tmpi = 0.0d0

           usr  = cshift(upxr,dim=xdir,shift=-ixq)
           usi  = cshift(upxi,dim=xdir,shift=-ixq)

           ustr = cshift(uptr,dim=xdir,shift=-ixq)
           usti = cshift(upti,dim=xdir,shift=-ixq)

           do ic = 1,nc
              do jc = 1,nc
                 do kc = 1,nc

                    tmpr(:,:,:,:,ic,jc) = tmpr(:,:,:,:,ic,jc)
     &                   + usr(:,:,:,:,kc,ic) * ustr(:,:,:,:,kc,jc)  + usi(:,:,:,:,kc,ic) * usti(:,:,:,:,kc,jc)

                    tmpi(:,:,:,:,ic,jc) = tmpi(:,:,:,:,ic,jc)  
     &                   + usr(:,:,:,:,kc,ic) * usti(:,:,:,:,kc,jc)  - usi(:,:,:,:,kc,ic) * ustr(:,:,:,:,kc,jc)

                 end do
              end do
           end do
c
c     go forwards from the negative x direction, forming staple
c
           sxnr = 0.0d0
           sxni = 0.0d0

           usr  = cshift(cshift(upxr,dim=that,shift= it),dim=xdir,shift=-ixq)
           usi  = cshift(cshift(upxi,dim=that,shift= it),dim=xdir,shift=-ixq)

           do ic = 1,nc
              do jc = 1,nc
                 do kc = 1,nc

                    sxnr(:,:,:,:,ic,jc) = sxnr(:,:,:,:,ic,jc)
     &                   + tmpr(:,:,:,:,ic,kc) * usr(:,:,:,:,kc,jc) - tmpi(:,:,:,:,ic,kc) * usi(:,:,:,:,kc,jc)

                    sxni(:,:,:,:,ic,jc) = sxni(:,:,:,:,ic,jc)
     &                   + tmpr(:,:,:,:,ic,kc) * usi(:,:,:,:,kc,jc) + tmpi(:,:,:,:,ic,kc) * usr(:,:,:,:,kc,jc)

                 end do
              end do
           end do
c
c-----------------------------------------------------------------------------------------
c-----------------------------------------------------------------------------------------
c-----------------------------------------------------------------------------------------
c
c  create staple in the positive y direction
c
           tmpr = 0.0d0
           tmpi = 0.0d0

           ustr = cshift(uptr,dim=ydir,shift= iyq)
           usti = cshift(upti,dim=ydir,shift= iyq)

           do ic = 1,nc
              do jc = 1,nc
                 do kc = 1,nc

                    tmpr(:,:,:,:,ic,jc) = tmpr(:,:,:,:,ic,jc)
     &                   + upyr(:,:,:,:,ic,kc) * ustr(:,:,:,:,kc,jc) - upyi(:,:,:,:,ic,kc) * usti(:,:,:,:,kc,jc)

                    tmpi(:,:,:,:,ic,jc) = tmpi(:,:,:,:,ic,jc)
     &                   + upyr(:,:,:,:,ic,kc) * usti(:,:,:,:,kc,jc) + upyi(:,:,:,:,ic,kc) * ustr(:,:,:,:,kc,jc)

                 end do
              end do
           end do
c
c     go backwards from the positive y direction, forming staple
c
           sypr = 0.0d0
           sypi = 0.0d0

           usr  = cshift(upyr,dim=that,shift= it)
           usi  = cshift(upyi,dim=that,shift= it)

           do ic = 1,nc
              do jc = 1,nc
                 do kc = 1,nc

                    sypr(:,:,:,:,ic,jc) = sypr(:,:,:,:,ic,jc)
     &                   + tmpr(:,:,:,:,ic,kc) * usr(:,:,:,:,jc,kc) + tmpi(:,:,:,:,ic,kc) * usi(:,:,:,:,jc,kc)

                    sypi(:,:,:,:,ic,jc) = sypi(:,:,:,:,ic,jc)
     &                   - tmpr(:,:,:,:,ic,kc) * usi(:,:,:,:,jc,kc) + tmpi(:,:,:,:,ic,kc) * usr(:,:,:,:,jc,kc)

                 end do
              end do
           end do
c
c---------------------------------------------------------------------------
c
c  create staple in the negative y direction
c
           tmpr = 0.0d0
           tmpi = 0.0d0

           usr  = cshift(upyr,dim=ydir,shift=-iyq)
           usi  = cshift(upyi,dim=ydir,shift=-iyq)

           ustr = cshift(uptr,dim=ydir,shift=-iyq)
           usti = cshift(upti,dim=ydir,shift=-iyq)

           do ic = 1,nc
              do jc = 1,nc
                 do kc = 1,nc

                    tmpr(:,:,:,:,ic,jc) = tmpr(:,:,:,:,ic,jc)
     &                   + usr(:,:,:,:,kc,ic) * ustr(:,:,:,:,kc,jc) + usi(:,:,:,:,kc,ic) * usti(:,:,:,:,kc,jc)

                    tmpi(:,:,:,:,ic,jc) = tmpi(:,:,:,:,ic,jc)
     &                   + usr(:,:,:,:,kc,ic) * usti(:,:,:,:,kc,jc) - usi(:,:,:,:,kc,ic) * ustr(:,:,:,:,kc,jc)

                 end do
              end do
           end do
c
c     go forwards from the negative y direction, forming staple
c
           synr = 0.0d0
           syni = 0.0d0

           usr  = cshift(cshift(upyr, dim=that,shift= it),dim=ydir,shift=-iyq)
           usi  = cshift(cshift(upyi, dim=that,shift= it),dim=ydir,shift=-iyq)

           do ic = 1,nc
              do jc = 1,nc
                 do kc = 1,nc

                    synr(:,:,:,:,ic,jc) = synr(:,:,:,:,ic,jc)
     &                   + tmpr(:,:,:,:,ic,kc) * usr(:,:,:,:,kc,jc) - tmpi(:,:,:,:,ic,kc) * usi(:,:,:,:,kc,jc)

                    syni(:,:,:,:,ic,jc) = syni(:,:,:,:,ic,jc)
     &                   + tmpr(:,:,:,:,ic,kc) * usi(:,:,:,:,kc,jc) + tmpi(:,:,:,:,ic,kc) * usr(:,:,:,:,kc,jc)

                 end do
              end do
           end do
