

pro PSD_Stitch_Plot, event
    widget_control,event.top,get_uvalue = pstate
    datasize = size(*((*pstate).data))
    if datasize(1) le 0 then begin
        void = dialog_message('The Data File is Empty!')
        return
    endif
    widget_control,(*pstate).PsdCenter, get_value=centerch
    widget_control,(*pstate).PsdLeft, get_value=leftch
    widget_control,(*pstate).PsdRight, get_value=rightch

    range=indgen(3)
    range(0) = leftch
    range(1) = centerch
    range(2) = rightch
    if range(1) lt 0 || range(1) gt 47 then begin
        void = dialog_message('The center of the PSD is out of [0 47], set to default value [0 23 47]')
        widget_control,(*pstate).PsdCenter,set_value = '23'
        widget_control,(*pstate).PsdLeft,set_value = '0'
        widget_control,(*pstate).PsdRight,set_value = '47'
        range(0) = 0
        range(1) = 23
        range(2) = 47
    endif
    if range(0) ge range(1) || range(0) ge range(2) || range(0) lt 0 || range(0) gt 47 then begin
        void = dialog_message('The left of the PSD is greater center, set default value to [0 23 47]')
        widget_control,(*pstate).PsdCenter,set_value = '23'
        widget_control,(*pstate).PsdLeft,set_value = '0'
        widget_control,(*pstate).PsdRight,set_value = '47'
        range(0) = 0
        range(1) = 23
        range(2) = 47
    endif
    if range(2) lt range(0) || range(2) lt range(1) || range(2) lt 0 || range(2) gt 47 then begin
        void = dialog_message('The left of the PSD is greater than the right, set default value to 0, 47')
        widget_control,(*pstate).PsdCenter,set_value = '23'
        widget_control,(*pstate).PsdLeft,set_value = '0'
        widget_control,(*pstate).PsdRight,set_value = '47'
        range(0) = 0
        range(1) = 23
        range(2) = 47
    endif

    Ch_eff=findgen(48)
    ;fn = strjoin([mydirectory,'PSD_Channeal_Eff.dat'])
    ;if fn ne '' then begin
;        openr, fun, fn, /get_lun
;        readf, fun, Ch_eff
;        free_lun, fun
;    endif
    Ch_eff=(*pstate).Ch_axis_Eff
    if (*pstate).colname eq 'A4' then begin
        data = (*(*pstate).data)
        var =(*(*pstate).vardata)
        sz = size(data)
        data_err = dindgen(sz(1), sz(2))
        data_err = sqrt(data)
        varname=*(*pstate).varname
        colname=(*pstate).colname
    var_pos=where(strlowcase(varname) eq strlowcase('A4'))
    var_pos=var_pos[0]
    ;Now open the A4 spacing calibration file and readout
    ch_position = fltarr(48)
    ;ch_space = dindgen(48)
    fn = (*pstate).PSD_A4_Spacing
    if (fn eq '') then a4_handler, event
    if fn eq '' then return 
    if fn ne '' then begin
      openr, fun, fn, /get_lun
      readf, fun, ch_position
      free_lun, fun
    endif
    ;Now calculate the spacing between every channel and the specified center channel
    ;for k = 0, 47 do begin
    ; ch_space(k) = ch_position(range(1))-ch_position(k)
    ;endfor
    ch_space=ch_position(range(1))-ch_position
    ch_space_extend=fltarr(49)
    ch_space_extend(0:47)=ch_space
    ch_space_extend(48)=ch_space_extend(47)
    ch_space=ch_space_extend
    ; After align the middle of PSD with the every A4 position and then build
    ; the left and right coordinate for every psd position
    ch_left = dblarr(49)
    ch_right =dblarr(48)
        A4_begin = var(var_pos, 0);-ch_space[range[0]]
        A4_end = var(var_pos, sz(2)-1);-ch_space[range[2]]
        if A4_begin eq A4_end then begin
            void = dialog_message('The step size of A4 is 0, could not stitch !')
            return
        endif
        output_width = 0.1
        output_npt = round(abs(A4_end-A4_begin)/output_width)
        output_data = dblarr(output_npt)
        data_norm=dblarr(output_npt)+1
        output_data_err = dblarr(output_npt)
        output_data_left = findgen(output_npt+1)
        output_data_right = findgen(output_npt)
        ;data_plus_count = findgen(output_npt)
        dis = dblarr(48, output_npt)
        frac = dblarr(48, output_npt)
        z_in = dblarr(48)
        dz_in = dblarr(48)
        output_data_left=min([A4_begin,A4_end])+output_data_left*output_width
        
        output_data_left(output_npt) = max([A4_begin,A4_end])
        output_data_right = output_data_left(1:output_npt)
        ;output_data_left=reverse(output_data_left)

        ;for i = 0, output_npt-1 do begin
        ;    data_plus_count(i) = 0.0
        ;endfor
        data_plus_count=dblarr(output_npt)
        mon_in=ch_left(0:47)+1.0
        for l=0,47 do begin
            if Ch_eff(l) lt 1E-2 then begin
            mon_in(l)=0
            endif
            endfor ; for
        output_tmp=data_plus_count
        dz_output_tmp=data_plus_count
        output_mon=data_plus_count
        dz_output_mon=data_plus_count
        ;print,'sz=',sz(2)
        widget_control, /hourglass
        for i = 0, sz(2)-1 do begin
            ;for j = 1, 47 do begin
            ;    ch_left(j) = ch_space(j) + 0.5*(ch_space(j-1) -ch_space(j))+var(0,i)
            ;endfor
            ch_left=ch_space+0.5*(-ch_space+shift(ch_space,1))+var(var_pos,i)
            ch_left(0) = ch_space(0) + 0.5*(ch_space(0) - ch_space(1)) +var(var_pos, i)
            ch_left(48) = ch_space(47) - 1.0*(ch_space(46)-ch_space(47)) + var(var_pos, i)
            ch_right = ch_left(1:48)
            z_in = data(*, i)*Ch_eff
            dz_in = data_err(*,i)*Ch_eff
            ;help,ch_left, z_in,dz_in,output_data_left
            ch_left=reverse(ch_left)
            z_in=reverse(z_in)
            dz_in=reverse(dz_in)
            ;output_data_left=reverse(output_data_left)
            drebin,ch_left[range[0]:range[2]+1],z_in[range[0]:range[2]],dz_in[range[0]:range[2]],output_data_left,output_tmp,dz_output_tmp,/histogram,/to_histogram,err=err,emsg=emsg
            print, emsg
            drebin,ch_left[range[0]:range[2]+1],mon_in[range[0]:range[2]],mon_in[range[0]:range[2]],output_data_left,output_mon,dz_output_mon,/histogram,/to_histogram,err=err,emsg=emsg
            output_data=output_data+output_tmp*output_mon
            output_data_err=output_data_err+dz_output_tmp^2*output_mon^2
            data_norm=data_norm+output_mon
            ;print, 'emsg2 ',emsg
            ;print,'output_mon ',n_elements(output_mon)
            ;print,'output_tmp ',n_elements(output_tmp)
            ;print,'ch_left',(ch_left)
            ;print,'data',n_elements(z_in)
            ;drebin_histo,x_in,z_in,dz_in,x_out,z_out,dz_out
            ;for j = range(0), range(2) do begin
            ;    for k = 0, output_npt-1 do begin
            ;        min_dis_righ = min([ch_left(j), output_data_right(k)])
            ;        max_dis_left = max([ch_right(j), output_data_left(k)])
            ;        dis(j, k) = max([0, (min_dis_righ-max_dis_left)])
            ;        frac(j, k) = dis(j, k)/abs(ch_left(j+1)-ch_left(j))
;
;                    output_data(k) = output_data(k)+z_in(j)*frac(j, k)
;                    output_data_err(k)= output_data_err(k) + (dz_in(j))^2*frac(j,k)
;                    if Ch_eff(j) gt 1E-2 then begin
;                        ;print,j,Ch_eff(j)
;                         data_plus_count(k) = data_plus_count(k) +frac(j, k)
;                    endif
;                endfor
;            endfor
        endfor; i
        output_data=output_data/data_norm
        output_data_err=sqrt(output_data_err)/data_norm
        ;for k = 0, output_npt-1 do begin
        ;    if data_plus_count(k) eq 0 then begin
        ;        print, 'Frac(', k,')=', 0
        ;        continue
        ;    endif
        ;    output_data(k) = output_data(k)/data_plus_count(k)
        ;    output_data_err(k)=sqrt(output_data_err(k))/data_plus_count(k)
        ;endfor
        xrdata = dindgen(3, output_npt)
        for k = 0, output_npt-1 do begin
            xrdata(0, k) = output_data_left(k)
            xrdata(1, k) = output_data(k)
            xrdata(2, k) = output_data_err(k)
        endfor
        *(*pstate).stitchdata = xrdata
        wset,(*pstate).winpix
        minx = min(xrdata(0,*), max=maxx)
        miny = min(xrdata(1,*), max=maxy)
        (*pstate).xdata = ptr_new(xrdata(0,*))
        (*pstate).ydata = ptr_new(xrdata(1,*))
        (*pstate).errdata = ptr_new(xrdata(2,*))
        (*pstate).plot_title = (*pstate).filename+': '+(*pstate).colname +' Stitch Plot' 
        (*pstate).xtitle = (*pstate).colname
        (*pstate).ytitle = 'Intensity ( '+ string(strtrim((*pstate).unit, 2))+' / '+(*pstate).unitname +' )'
        (*pstate).plottype =0
    PSD_Plot_graphs, pstate
        
        wset,(*pstate).winvis
    device, copy=[0,0,700,700,0,0,(*pstate).winPix] ; 
    endif; for A4