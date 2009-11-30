      sbxlr = cshift( bxlr,dim=that,shift=it)
      sbxli = cshift( bxli,dim=that,shift=it)

      sbllr = cshift( bllr,dim=that,shift=it)
      sblli = cshift( blli,dim=that,shift=it)

      sbrlr = cshift( brlr,dim=that,shift=it)
      sbrli = cshift( brli,dim=that,shift=it)

      Res = 0.0d0
      do ic = 1,36
         do jc = 1,36
            do kc = 1,36

               tmp1r(:,:,:,:) =
     &              bxlr(:,:,:,:,la(ic),lb(jc)) * bllr(:,:,:,:,lb(kc),la(jc)) * brlr(:,:,:,:,la(kc),lb(ic))
     &            + bxlr(:,:,:,:,la(ic),lb(jc)) * blli(:,:,:,:,lb(kc),la(jc)) * brli(:,:,:,:,la(kc),lb(ic))
     &            - bxli(:,:,:,:,la(ic),lb(jc)) * bllr(:,:,:,:,lb(kc),la(jc)) * brli(:,:,:,:,la(kc),lb(ic))
     &            + bxli(:,:,:,:,la(ic),lb(jc)) * blli(:,:,:,:,lb(kc),la(jc)) * brlr(:,:,:,:,la(kc),lb(ic))

               tmp1i(:,:,:,:) =
     &            + bxli(:,:,:,:,la(ic),lb(jc)) * blli(:,:,:,:,lb(kc),la(jc)) * brli(:,:,:,:,la(kc),lb(ic))
     &            + bxlr(:,:,:,:,la(ic),lb(jc)) * bllr(:,:,:,:,lb(kc),la(jc)) * brli(:,:,:,:,la(kc),lb(ic))
     &            - bxlr(:,:,:,:,la(ic),lb(jc)) * blli(:,:,:,:,lb(kc),la(jc)) * brlr(:,:,:,:,la(kc),lb(ic))
     &            + bxli(:,:,:,:,la(ic),lb(jc)) * bllr(:,:,:,:,lb(kc),la(jc)) * brlr(:,:,:,:,la(kc),lb(ic))

               tmp2r(:,:,:,:) =
     &              tl1r(:,:,:,:,lc(jc),lcp(jc)) * tl2r(:,:,:,:,lc(kc),lcp(kc)) * tlr(:,:,:,:,lc(ic),lcp(ic))
     &            - tl1r(:,:,:,:,lc(jc),lcp(jc)) * tl2i(:,:,:,:,lc(kc),lcp(kc)) * tli(:,:,:,:,lc(ic),lcp(ic))
     &            - tl1i(:,:,:,:,lc(jc),lcp(jc)) * tl2r(:,:,:,:,lc(kc),lcp(kc)) * tli(:,:,:,:,lc(ic),lcp(ic))
     &            - tl1i(:,:,:,:,lc(jc),lcp(jc)) * tl2i(:,:,:,:,lc(kc),lcp(kc)) * tlr(:,:,:,:,lc(ic),lcp(ic))

               tmp2i(:,:,:,:) =
     &            - tl1i(:,:,:,:,lc(jc),lcp(jc)) * tl2i(:,:,:,:,lc(kc),lcp(kc)) * tli(:,:,:,:,lc(ic),lcp(ic))
     &            + tl1r(:,:,:,:,lc(jc),lcp(jc)) * tl2r(:,:,:,:,lc(kc),lcp(kc)) * tli(:,:,:,:,lc(ic),lcp(ic))
     &            + tl1r(:,:,:,:,lc(jc),lcp(jc)) * tl2i(:,:,:,:,lc(kc),lcp(kc)) * tlr(:,:,:,:,lc(ic),lcp(ic))
     &            + tl1i(:,:,:,:,lc(jc),lcp(jc)) * tl2r(:,:,:,:,lc(kc),lcp(kc)) * tlr(:,:,:,:,lc(ic),lcp(ic))

               tmp3r(:,:,:,:) =
     &              sbxlr(:,:,:,:,lap(ic),lbp(jc)) * sbllr(:,:,:,:,lbp(kc),lap(jc)) * sbrlr(:,:,:,:,lap(kc),lbp(ic))
     &            + sbxlr(:,:,:,:,lap(ic),lbp(jc)) * sblli(:,:,:,:,lbp(kc),lap(jc)) * sbrli(:,:,:,:,lap(kc),lbp(ic))
     &            - sbxli(:,:,:,:,lap(ic),lbp(jc)) * sbllr(:,:,:,:,lbp(kc),lap(jc)) * sbrli(:,:,:,:,lap(kc),lbp(ic))
     &            + sbxli(:,:,:,:,lap(ic),lbp(jc)) * sblli(:,:,:,:,lbp(kc),lap(jc)) * sbrlr(:,:,:,:,lap(kc),lbp(ic))

               tmp3i(:,:,:,:) =
     &            + sbxli(:,:,:,:,lap(ic),lbp(jc)) * sblli(:,:,:,:,lbp(kc),lap(jc)) * sbrli(:,:,:,:,lap(kc),lbp(ic))
     &            + sbxlr(:,:,:,:,lap(ic),lbp(jc)) * sbllr(:,:,:,:,lbp(kc),lap(jc)) * sbrli(:,:,:,:,lap(kc),lbp(ic))
     &            - sbxlr(:,:,:,:,lap(ic),lbp(jc)) * sblli(:,:,:,:,lbp(kc),lap(jc)) * sbrlr(:,:,:,:,lap(kc),lbp(ic))
     &            + sbxli(:,:,:,:,lap(ic),lbp(jc)) * sbllr(:,:,:,:,lbp(kc),lap(jc)) * sbrlr(:,:,:,:,lap(kc),lbp(ic))

               Res(:,:,:,:) = Res(:,:,:,:) + fper(ic) * fper(jc) * fper(kc) * (
     &              tmp1r(:,:,:,:) * tmp2r(:,:,:,:) * tmp3r(:,:,:,:)
     &            - tmp1r(:,:,:,:) * tmp2i(:,:,:,:) * tmp3i(:,:,:,:)
     &            - tmp1i(:,:,:,:) * tmp2r(:,:,:,:) * tmp3i(:,:,:,:)
     &            - tmp1i(:,:,:,:) * tmp2i(:,:,:,:) * tmp3r(:,:,:,:) )

            end do
         end do
      end do
