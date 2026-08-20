# get datasets from repo
# try simple load, download if fail
rdataload <- function(dataurl = NULL) {
  x <- gsub('\\.RData', '', basename(dataurl))

  # try simple load
  ld <- try(load(url(dataurl)), silent = T)

  # return x if load worked
  if (!inherits(ld, 'try-error')) {
    out <- get(x)
  }

  # download x if load failed
  if (inherits(ld, 'try-error')) {
    fl <- paste(tempdir(), basename(dataurl), sep = '/')
    download.file(flurl, destfile = fl, quiet = T)
    load(file = fl)
    out <- get(x)
    suppressMessages(file.remove(fl))
  }

  return(out)
}

# function for formatting p-values in tables
p_ast <- function(x) {
  sig_cats <- c('**', '*', '')
  sig_vals <- c(-Inf, 0.005, 0.05, Inf)

  out <- cut(x, breaks = sig_vals, labels = sig_cats, right = FALSE)
  out <- as.character(out)

  return(out)
}

# p-value text formatting for GAM summary table
p_txt <- function(x, addp = TRUE) {
  if (x < 0.001) {
    out <- 'p < 0.001'
  } else {
    out <- paste('p =', sprintf('%.3f', round(x, 3)))
  }

  if (!addp) {
    out <- gsub('^p\\s|\\=\\s', '', out)
  }

  return(out)
}

# model rsq, same as summary(mod)$dev.expl
rsq_fun <- function(mod) {
  # get complete cases
  toeval <- data.frame(resid = mod$residuals, obs = mod$y)
  toeval <- na.omit(toeval)

  ssr <- sum(toeval$resid^2)
  sst <- sum((toeval$obs - mean(toeval$obs))^2)

  out <- 1 - (ssr / sst)

  return(out)
}
# get salinity prediction grids
fits_fun <- function(modin, toprd, salgrd, type = 'link') {
  tibble(
    res = as.numeric(predict(modin, newdata = toprd, type = type))
  ) |>
    bind_cols(toprd) |>
    mutate(
      year = lubridate::year(date),
      month = lubridate::month(date),
      day = lubridate::day(date),
      sal = factor(
        sal,
        levels = salgrd,
        labels = paste0('X', seq(1:length(salgrd)))
      )
    ) |>
    select(-dec_time, -doy) |>
    pivot_wider(names_from = sal, values_from = res)
}

# get gam predictions
pred_fun <- function(datin, modin) {
  moddat <- datin
  moddat$fit <- predict(modin, type = 'link', newdata = moddat)
  moddat$btfit <- predict(modin, type = 'response', newdata = moddat)
  moddat$fitmd <- predict(
    modin,
    type = 'link',
    newdata = moddat %>% mutate(tn_load = mean(tn_load))
  )
  moddat$btfitmd <- predict(
    modin,
    type = 'response',
    newdata = moddat %>% mutate(tn_load = mean(tn_load))
  )

  salgrd <- moddat |>
    pull(sal) |>
    range(na.rm = TRUE)
  salgrd <- seq(salgrd[1], salgrd[2], length.out = 10)

  # make prediction grids
  toprd <- moddat |>
    select(date, dec_time, doy, tn_load) |>
    crossing(
      sal = salgrd
    )

  toprdmd <- moddat |>
    select(date, dec_time, doy, tn_load) |>
    mutate(tn_load = mean(tn_load)) |>
    crossing(
      sal = salgrd
    )

  # get link predictions, wide format
  fits <- fits_fun(modin, toprd, salgrd, type = 'link')
  fitsmd <- fits_fun(modin, toprdmd, salgrd, type = 'link')

  # get response predictions, wide format
  btfits <- fits_fun(modin, toprd, salgrd, type = 'response')
  btfitsmd <- fits_fun(modin, toprdmd, salgrd, type = 'response')

  # normalized results
  norm <- norm_fun(moddat, fits, btfits, salgrd)
  normmd <- norm_fun(moddat, fitsmd, btfitsmd, salgrd)[, c('norm', 'btnorm')] |>
    rename(normmd = norm, btnormmd = btnorm)

  moddat <- bind_cols(norm, normmd)

  out <- structure(
    .Data = moddat,
    class = c('data.frame', 'tibble'),
    fits = fits,
    btfits = btfits,
    fitsmd = fitsmd,
    btfitsmd = btfitsmd,
    salgrd = salgrd
  )

  return(out)
}

# get annual summaries of predictions
annsum_fun <- function(prds) {
  out <- data.frame(prds) |>
    mutate(
      yr = lubridate::year(date)
    ) |>
    summarise(
      chla = mean(chla, na.rm = T),
      btfit = mean(btfit, na.rm = T),
      btnorm = mean(btnorm, na.rm = T),
      btfitmd = mean(btfitmd, na.rm = T),
      btnormmd = mean(btnormmd, na.rm = T),
      .by = c(yr)
    )

  return(out)
}

# gridplot for gam
# ... additional arguments to pass to gridpred_fun
grid_plo <- function(
  prds,
  month = c(1:12),
  years = NULL,
  col_vec = NULL,
  col_lim = NULL,
  allsal = FALSE,
  sal_fac = 3,
  yr_fac = 3,
  thresh = NULL,
  ncol = NULL,
  grids = FALSE,
  pretty = TRUE,
  ldmod = 'btfits',
  ...
) {
  # convert month vector to those present in data
  allmo <- FALSE
  if ('all' %in% month) {
    allmo <- TRUE
    month <- c(1:12)
  }

  ldmod <- match.arg(
    ldmod,
    choices = c('btfits', 'btfitsmd')
  )

  # format predictions as wide
  to_plo <- attr(prds, ldmod)

  # model data with date column
  moddat <- data.frame(prds) |>
    filter(
      lubridate::year(date) >= min(to_plo$year) &
        lubridate::year(date) <= max(to_plo$year)
    )

  # salinity grid values
  salgrd <- attr(prds, 'salgrd')

  # convert month vector to those present in data
  month <- month[month %in% to_plo$month]
  if (length(month) == 0) {
    stop('No observable data for the chosen month')
  }

  # axis labels
  ylabel <- 'Chlorophyll-a (µg/L)'
  xlabel <- 'Salinity (ppth)'

  # subset years to plot
  if (!is.null(years)) {
    if (length(years) != 2) {
      stop('years argument must have two values for first and last')
    }
  } else {
    years <- range(to_plo$year)
  }

  years <- seq(years[1], years[2])
  to_plo <- to_plo[to_plo$year %in% years, ]

  if (nrow(to_plo) == 0) {
    stop('No data to plot for the date range')
  }

  # reshape data frame
  to_plo <- to_plo[to_plo$month %in% month, , drop = FALSE]
  names(to_plo)[grep('^X', names(to_plo))] <- paste('sal', salgrd)
  to_plo <- tidyr::gather(to_plo, 'sal', 'res', 6:ncol(to_plo)) |>
    mutate(sal = as.numeric(gsub('^sal ', '', sal))) |>
    select(-date, -day) |>
    summarize(
      res = mean(res, na.rm = TRUE),
      .by = c(year, month, sal)
    )

  ## use linear interpolation to make a smoother plot
  if (!allmo) {
    # these are factors by which salinity and years are multiplied for interpolation
    sal_fac <- length(salgrd) * sal_fac
    sal_fac <- seq(min(salgrd), max(salgrd), length.out = sal_fac)
    yr_fac <- length(unique(to_plo$year)) * yr_fac
    yr_fac <- seq(min(to_plo$year), max(to_plo$year), length.out = yr_fac)

    # separately by month
    to_plo <- split(to_plo, to_plo$month)

    to_plo <- lapply(to_plo, function(x) {
      # interp across salinity first
      interped <- lapply(
        split(x, x$year),
        function(y) {
          out <- approx(y$sal, y$res, xout = sal_fac, rule = 2)
          out <- data.frame(year = unique(y$year), month = unique(y$month), out)
          return(out)
        }
      )
      interped <- do.call('rbind', interped)
      names(interped) <- c('year', 'month', 'sal', 'res')

      # interp across years
      interped <- lapply(
        split(interped, interped$sal),
        function(y) {
          out <- approx(y$year, y$res, xout = yr_fac, rule = 2)
          out <- data.frame(
            year = out$x,
            month = unique(y$month),
            sal = unique(y$sal),
            res = out$y
          )
          return(out)
        }
      )
      interped <- do.call('rbind', interped)
      names(interped) <- c('year', 'month', 'sal', 'res')

      return(interped)
    })

    to_plo <- do.call('rbind', to_plo)
    row.names(to_plo) <- 1:nrow(to_plo)
  }

  ## use linear interpolation to make a smoother plot, handle differently if allmo
  if (allmo) {
    # format to_plo for interp (wide, as matrix)
    to_interp <- to_plo
    to_interp$date <- with(
      to_interp,
      as.Date(paste(year, month, '1', sep = '-'))
    )
    to_interp <- ungroup(to_interp) |>
      select(date, sal, res) |>
      tidyr::spread(sal, res)

    # values to pass to interp
    dts <- lubridate::decimal_date(to_interp$date)
    fit_grd <- select(to_interp, -date)
    sal_fac <- length(salgrd) * sal_fac
    sal_fac <- seq(min(salgrd), max(salgrd), length.out = sal_fac)
    yr_fac <- seq(min(dts), max(dts), length.out = length(dts) * yr_fac)
    to_norm <- expand.grid(yr_fac, sal_fac)

    # bilinear interpolation of fit grid with data to average for norms
    norms <- fields::interp.surface(
      obj = list(
        y = salgrd,
        x = dts,
        z = data.frame(fit_grd)
      ),
      loc = to_norm
    )

    to_plo <- data.frame(to_norm, norms)
    names(to_plo) <- c('year', 'sal', 'res')
  }

  # constrain plots to salinity limits for the selected month
  if (!allsal & !allmo) {
    #min, max salinity values to plot
    lim_vals <- moddat |>
      mutate(
        month = lubridate::month(date),
        yr = lubridate::year(date)
      ) |>
      filter(yr %in% years) |>
      summarize(
        Low = quantile(sal, 0, na.rm = TRUE),
        High = quantile(sal, 1, na.rm = TRUE),
        .by = month
      )

    # month sal ranges for plot
    lim_vals <- lim_vals[lim_vals$month %in% month, ]

    # merge limts with months
    to_plo <- left_join(to_plo, lim_vals, by = 'month')

    # reduce data
    sel_vec <- with(
      to_plo,
      sal >= Low &
        sal <= High
    )
    to_plo <- to_plo[sel_vec, !names(to_plo) %in% c('Low', 'High')]
    to_plo <- arrange(to_plo, year, month)
  }

  # contstrain all data by quantiles if not separated by month
  if (!allsal & allmo) {
    quants <- moddat |>
      pull(sal) |>
      quantile(c(0.05, 0.95), na.rm = TRUE)
    to_plo <- to_plo[with(to_plo, sal >= quants[1] & sal <= quants[2]), ]
  }

  # change month vector of not plotting all months in same plot
  if (!allmo) {
    # months labels as text
    mo_lab <- data.frame(
      num = seq(1:12),
      txt = c(
        'January',
        'February',
        'March',
        'April',
        'May',
        'June',
        'July',
        'August',
        'September',
        'October',
        'November',
        'December'
      )
    )
    mo_lab <- mo_lab[mo_lab$num %in% month, ]
    to_plo$month <- factor(
      to_plo$month,
      levels = mo_lab$num,
      labels = mo_lab$txt
    )
  }

  # make plot
  p <- ggplot(to_plo, aes(x = year, y = sal, fill = res)) +
    geom_tile(data = subset(to_plo, !is.na(to_plo$res)), aes(fill = res)) +
    geom_tile(
      data = subset(to_plo, is.na(to_plo$res)),
      fill = 'black',
      alpha = 0
    )

  if (thresh) {
    p <- p +
      geom_contour(
        aes(z = res),
        breaks = thresh,
        color = "grey",
        linetype = "solid",
        linewidth = 1
      )
  }

  if (!allmo) {
    p <- p + facet_wrap(~month, ncol = ncol)
  }

  # return bare bones if FALSE
  if (!pretty) {
    return(p)
  }

  # get colors
  cols <- cols <- RColorBrewer::brewer.pal(11, 'Spectral')

  p <- p +
    theme_minimal() +
    theme(
      legend.position = 'top',
      axis.title.x = element_blank(),
      plot.title = element_text(hjust = 0.5)
    ) +
    # scale_x_continuous(expand = c(0.05, 0.05)) +
    scale_y_continuous(xlabel) + #, expand = c(0.05,0.05)) +
    scale_fill_gradientn(ylabel, colours = rev(cols), limits = col_lim) +
    guides(fill = guide_colourbar(barwidth = 10))

  # add grid lines
  if (!grids) {
    p <- p +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
  }

  return(p)
}

# normalize predictions
norm_fun <- function(dat_in, fits, btfits, salgrd) {
  num_obs <- nrow(dat_in)

  # prep interp grids by adding month, year columns
  dts <- fits$date
  fits <- select(fits, -year, -month, -day, -date, -tn_load)
  btfits <- btfits |>
    select(-year, -date, -month, -day, -tn_load)

  # sal values occuring by month, used for interpolation
  sal_mon <- data.frame(dat_in) |>
    select(date, sal) |>
    mutate(
      yr = lubridate::year(date),
      mo = lubridate::month(date)
    ) |>
    select(-date) |>
    summarize(
      sal = mean(sal, na.rm = TRUE),
      .by = c(yr, mo)
    ) |>
    tidyr::spread(yr, sal)

  # values to interpolate for normalization
  to_norm <- data.frame(dat_in) |>
    select(date) |>
    mutate(mo = lubridate::month(date)) |>
    left_join(sal_mon, by = 'mo') |>
    tidyr::gather('yr', 'sal', -date, -mo) |>
    arrange(date) |>
    select(date, sal) |>
    na.omit()

  # bilinear interpolatoin of fit grid with data to average for norms
  norms <- fields::interp.surface(
    obj = list(
      y = salgrd,
      x = dts,
      z = data.frame(fits)
    ),
    loc = to_norm
  )

  # append to normalization data, then average for unique dates
  # averaging happens only if more than 80% of the predictions for a date are filled
  norms <- data.frame(date = to_norm$date, norms) |>
    summarise(
      norm = ifelse(
        sum(is.na(norms)) / length(norms) > 0.2,
        NA,
        mean(norms, na.rm = TRUE)
      ),
      .by = date
    )

  # bilinear interpolatoin of fit grid with data to average for norms
  btnorms <- fields::interp.surface(
    obj = list(
      y = salgrd,
      x = dts,
      z = data.frame(btfits)
    ),
    loc = to_norm
  )

  # append to normalization data, then average for unique dates
  # averaging happens only if more than 80% of the predictions for a date are filled
  btnorms <- data.frame(date = to_norm$date, btnorms) |>
    summarise(
      btnorm = ifelse(
        sum(is.na(btnorms)) / length(btnorms) > 0.2,
        NA,
        mean(btnorms, na.rm = TRUE)
      ),
      .by = date
    )

  # append to dat_in object
  dat_in <- dplyr::left_join(dat_in, norms, by = 'date')
  dat_in <- dplyr::left_join(dat_in, btnorms, by = 'date')

  # exit function
  return(dat_in)
}

# add lags to load data, then get cumulative loads by lag
lddatlag_fun <- function(lddat) {
  out <- lddat |>
    summarise(
      tn_load = sum(tn_load),
      .by = c(bay_segment, date)
    ) |>
    mutate(
      lag1 = dplyr::lag(tn_load, n = 1),
      lag2 = dplyr::lag(tn_load, n = 2),
      lag3 = dplyr::lag(tn_load, n = 3),
      .by = c(bay_segment)
    ) |>
    fill(
      lag1,
      lag2,
      lag3,
      .by = c(bay_segment),
      .direction = 'up'
    ) |>
    mutate(
      tn_load0 = tn_load,
      tn_load1 = tn_load + lag1,
      tn_load2 = tn_load + lag1 + lag2,
      tn_load3 = tn_load + lag1 + lag2 + lag3
    ) |>
    select(-tn_load, -lag1, -lag2, -lag3)

  return(out)
}

# prep data for modeling with load lags
# used for addl GAMs evaluating residuals, lags, and flexible k
modprep_fun <- function(wqdat, lddat) {
  lddatlag <- lddatlag_fun(lddat)

  wqdat <- wqdat |>
    mutate(salraw = sal) |>
    rename(
      tn_loadraw = tn_load
    ) |>
    left_join(lddatlag, by = c('bay_segment', 'date'))

  out <- wqdat |>
    filter(dec_time >= 1985) |>
    pivot_longer(
      cols = ends_with(c('load0', 'load1', 'load2', 'load3')),
      names_to = 'tn_load_lag',
      values_to = 'tn_load'
    ) |>
    group_nest(bay_segment, tn_load_lag)

  return(out)
}

# get load scalars for forecasting
# loads are scaled by ldfac based on monthly proportion in each year from historical record
# then summed by lags
ldscale_fun <- function(lddat, ldfac = c(0.5, 1, 2), yrs, bay_segment) {
  lds <- tibble(
    bay_segment = c('OTB', 'HB', 'MTB', 'LTB'),
    thresh = c(486, 1451, 799, 349)
  )

  perc <- lddat |>
    filter(bay_segment %in% !!bay_segment) |>
    select(-bay_segment) |>
    filter(yr %in% yrs) |>
    mutate(
      perc = tn_load / sum(tn_load),
      thresh = lds$thresh[lds$bay_segment == !!bay_segment],
      .by = yr
    )

  actload <- perc |>
    arrange(yr, mo) |>
    mutate(
      ldfac = 'Actual Load',
      tn_loadann = sum(tn_load),
      .by = yr
    )

  out <- perc |>
    crossing(
      ldfac = ldfac
    ) |>
    mutate(
      tn_loadann = thresh * ldfac,
      tn_load = perc * tn_loadann
    ) |>
    arrange(ldfac, yr, mo) |>
    filter(yr %in% yrs) |>
    mutate(
      ldfac = paste0(
        'TMDL Load Factor: ',
        ldfac
      )
    ) |>
    bind_rows(actload) |>
    mutate(
      ldfac = factor(
        ldfac,
        levels = sort(unique(ldfac))
      )
    )

  return(out)
}

# prep simulation data for salinity/load change scenarios
simprp_fun <- function(wqdat, lddat, yrs = c(2017:2021)) {
  # check if yrs arg sequential
  if (any(diff(yrs) > 1)) {
    stop('yrs argument must be sequential')
  }

  # calculate annual salinity change rates
  salyrrates <- wqdat |>
    dplyr::select(date, sal, bay_segment) |>
    mutate(
      yr = year(date)
    ) |>
    summarise(
      sal = mean(sal, na.rm = T),
      .by = c(yr, bay_segment)
    ) |>
    group_nest(bay_segment) |>
    mutate(
      salmod = purrr::map(
        data,
        ~ lm(sal ~ yr, data = .x, na.action = 'na.exclude')
      ),
      slopeyr = purrr::map_dbl(salmod, ~ coef(.x)[2])
    ) |>
    select(-salmod, -data) |>
    unnest(c('slopeyr')) |>
    crossing(
      yrs = c(1:50)
    ) |>
    mutate(
      futureyr = 2024 + yrs,
      salforeyr = slopeyr * yrs,
      .by = bay_segment
    )

  # get salinity avgs by yr, mo for test period
  saltstyrmo <- wqdat |>
    dplyr::select(date, sal, bay_segment) |>
    mutate(
      yr = year(date),
      mo = month(date)
    ) |>
    filter(yr %in% yrs) |>
    summarise(
      sal = mean(sal, na.rm = T),
      .by = c(yr, mo, bay_segment)
    )

  # add salinity yr, mo avg for test period to forecasts
  salrates <- salyrrates |>
    left_join(
      saltstyrmo,
      by = c('bay_segment'),
      relationship = 'many-to-many'
    ) |>
    mutate(
      salforeyrmo = sal + salforeyr,
      .by = bay_segment
    ) |>
    mutate(
      future_mo = mo,
      date = make_date(year = yr, month = mo, day = 1),
      doy = yday(date),
      dec_time = decimal_date(date)
    )

  # # add doy correction to salinity - skipped
  # salrates <- salrates |>
  #   group_nest(bay_segment, yrs) |>
  #   mutate(
  #     data = purrr::map(data, function(x) {
  #       salmod <- gam(
  #         salforeyrmo ~ s(doy, bs = 'cc'),
  #         data = x,
  #         method = 'REML',
  #         knots = list(doy = c(0, 366))
  #       )
  #       x$salforeyrmo <- residuals(salmod)
  #       return(x)
  #     })
  #   ) |>
  #   unnest('data')

  # get load scaling specific to yrs
  ldscale <- tibble(
    bay_segment = c('OTB', 'HB', 'MTB', 'LTB')
  ) |>
    mutate(
      ldscale = purrr::map(
        bay_segment,
        ~ ldscale_fun(
          lddat,
          ldfac = c(0.5, 1, 2),
          yrs = yrs,
          bay_segment = .x
        )
      )
    ) |>
    unnest('ldscale')

  # # add lagged loadings - skip
  # ldlag <- ldscale |>
  #   group_nest(ldfac) |>
  #   mutate(
  #     data = purrr::map(data, function(x) {
  #       lddatlag_fun(x)
  #     })
  #   ) |>
  #   unnest('data') |>
  #   select(ldfac, bay_segment, date, tn_load3)

  # # get year range and correct lag - skip
  # ldscale <- ldscale |>
  #   left_join(
  #     ldlag,
  #     by = c('bay_segment', 'date', 'ldfac'),
  #     relationship = 'many-to-many'
  #   ) |>
  #   filter(yr %in% yrs) |>
  #   mutate(
  #     tn_load = tn_load3,
  #     .keep = 'unused'
  #   )

  # join salinity projections with load scenarios
  # correct load for salinity and doy - skipped
  out <- salrates |>
    arrange(bay_segment, yrs, date) |>
    left_join(
      ldscale,
      by = c('bay_segment', 'date', 'yr', 'mo'),
      relationship = 'many-to-many'
    ) |>
    select(-sal) |>
    rename(sal = salforeyrmo) |>
    # group_nest(bay_segment, yrs, ldfac) |>
    # mutate(
    #   data = purrr::map(data, function(x) {
    #     tnmod <- gam(
    #       tn_load ~ s(doy, bs = 'cc') + s(sal),
    #       data = x,
    #       method = 'REML',
    #       knots = list(doy = c(0, 366))
    #     )
    #     x$tn_load <- residuals(tnmod)
    #     x
    #   })
    # ) |>
    # unnest(data) |>
    arrange(bay_segment, ldfac, yrs, date)

  return(out)
}

# simulation predictions for scenarios
simprd_fun <- function(modin, simdat, nsims = 100, chunk = T, all = F) {
  # modin <- modin |>
  #   filter(tn_load_lag == 'tn_load3') |>
  #   select(-tn_load_lag)

  trgs <- tbeptools::targets |>
    filter(bay_segment %in% modin$bay_segment) |>
    dplyr::select(bay_segment, thresh = chla_thresh) |>
    mutate(
      bay_segment = factor(bay_segment, levels = c('OTB', 'HB', 'MTB', 'LTB'))
    ) |>
    tibble()

  tojn <- simdat |>
    group_nest(bay_segment) |>
    rename(tst = data) |>
    left_join(trgs, by = c('bay_segment'), relationship = 'one-to-one')

  out <- modin |>
    left_join(tojn, by = c('bay_segment'), relationship = 'one-to-one')

  if (chunk) {
    out <- out |>
      mutate(
        simsyr = purrr::pmap(
          list(mod, tst, bay_segment),
          function(mod, tst, bay_segment) {
            out <- NULL
            simi <- nsims / 10
            reps <- rep(simi, 10)
            # chunking needs to be done to avoid memory issues
            for (i in seq_along(reps)) {
              cat(bay_segment, i, 'of', length(reps), '\n')
              simout <- simulate(mod, data = tst, nsim = simi) |>
                bind_cols(tst) |>
                pivot_longer(
                  cols = starts_with('sim'),
                  names_to = 'sim',
                  values_to = 'chla_sim'
                ) |>
                mutate(
                  sim = i * as.integer(str_remove(sim, 'sim_'))
                ) |>
                arrange(sim, date) |>
                summarise(
                  chla_sim = mean(chla_sim, na.rm = T),
                  .by = c(yrs, yr, ldfac, sim)
                )

              out <- bind_rows(out, simout)
            }

            return(out)
          }
        )
      )
  }

  if (!chunk) {
    out <- out |>
      mutate(
        sims = purrr::pmap(list(mod, tst), function(mod, tst) {
          simulate(mod, data = tst, nsim = nsims) |>
            bind_cols(tst) |>
            pivot_longer(
              cols = starts_with('sim'),
              names_to = 'sim',
              values_to = 'chla_sim'
            ) |>
            mutate(
              sim = as.integer(str_remove(sim, 'sim_'))
            ) |>
            arrange(sim, date)
        }),
        simsyr = purrr::map(sims, function(sims) {
          sims |>
            summarise(
              chla_sim = mean(chla_sim, na.rm = T),
              .by = c(yrs, yr, ldfac, sim)
            )
        })
      ) |>
      select(-sims)
  }

  out <- out |>
    mutate(
      exceedsyr = purrr::pmap(list(simsyr, thresh), function(simsyr, thresh) {
        simsyr |>
          mutate(
            exceeds = chla_sim > thresh
          ) |>
          summarise(
            percexceeds = mean(exceeds, na.rm = T) * 100,
            .by = c(yrs, sim, ldfac)
          )
      }),
      exceedssum = purrr::map(exceedsyr, function(exceedsyr) {
        exceedsyr |>
          summarise(
            avexceeds = mean(percexceeds, na.rm = T),
            sdexceeds = sd(percexceeds, na.rm = T),
            .by = c(yrs, ldfac)
          )
      })
    )

  if (!all) {
    out <- out |>
      select(-data, -mod, -prds, -annsum, -simsyr)
  }

  return(out)
}
