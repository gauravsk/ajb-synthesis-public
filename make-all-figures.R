# R code for generating figures in 
# "Quantifying soil microbial effects on plant species coexistence: 
# a conceptual synthesis", submitted to Am. J. Botany
# by Gaurav Kandlikar

requirements <- c("tidyverse", "patchwork", "deSolve", "rootSolve", "latex2exp")
# Load the required packages, or if a package is not available, install from CRAN
lapply(requirements, function(x) {
  if(!require(x, character.only = T)){
    install.packages(x)
  }
})


# Define a function used for making the PSF framework
# schematic, given a set of parameter values
make_params_plot <- function(params, scale = 1.5) {
  
  color_func <- function(x) {
    ifelse(x < 0, "#EE6677", "#4477AA")
  }
  df <- data.frame(x = c(0,0,1,1),
                   y = c(0,1,0,1),
                   type = c("M", "P", "M", "P"))
  
  params_plot <- 
    ggplot(df) +
    geom_label(x = 0, y = 1.1, size = 3.25, label = "Plant 1", color = "white", fill = "#9970ab", fontface = "bold") + 
    geom_label(x = 1, y = 1.1, size = 3.25, label = "Plant 2", color = "white", fill = "#5aae61", fontface = "bold") +
    geom_label(x = 0, y = -0.15, size = 3.25, label = "Soil\nmicrobes A", label.size = 0) + 
    geom_label(x = 1, y = -0.15, size = 3.25, label = "Soil\nmicrobes B", label.size = 0) +
    geom_label(x = 0.45, y = -0.4, size = 3.00, label.size = 0.05,
               label = TeX(paste0("$I_S$ = (",
                                  params["m1A"], " + ",
                                  params["m2B"], ") - (",
                                  params["m1B"], " + ", 
                                  params["m2A"], ") = ",
                                  params["m1A"] + params["m2B"] - params["m1B"] - params["m2A"]))) + 
    geom_segment(aes(x = 0, xend = 0, y = 0.1, yend = 0.9),
                 arrow = arrow(length = unit(0.03, "npc")),
                 linewidth = abs(params["m1A"])*scale, 
                 color = alpha(color_func(params["m1A"]), 1)) +
    geom_segment(aes(x = 0.05, xend = 0.95, y = 0.1, yend = 0.9),
                 arrow = arrow(length = unit(0.03, "npc")),
                 linewidth = abs(params["m2A"])*scale, 
                 color = alpha(color_func(params["m1B"]),1)) +
    geom_segment(aes(x = 0.95, xend = 0.05, y = 0.1, yend = 0.9),
                 arrow = arrow(length = unit(0.03, "npc")),
                 linewidth = abs(params["m1B"])*scale, 
                 color = alpha(color_func(params["m2A"]), 1)) +
    geom_segment(aes(x = 1, xend = 1, y = 0.1, yend = 0.9),
                 arrow = arrow(length = unit(0.03, "npc"),),
                 linewidth = abs(params["m2B"])*scale, 
                 color = alpha(color_func(params["m2B"]), 1)) +
    
    # Plant cultivation of microbes
    geom_segment(aes(x = -0.25, xend = -0.25, y = 0.9, yend = 0.1), linewidth = 0.15, linetype = 1,
                 arrow = arrow(length = unit(0.03, "npc"))) +
    geom_segment(aes(x = 1.25, xend = 1.25, y = 0.9, yend = 0.1), linewidth = 0.15, linetype = 1,
                 arrow = arrow(length = unit(0.03, "npc"))) +
    
    annotate("text", x = 0, y = 0.5, 
             label = TeX(paste0("m1A = ", params["m1A"])), 
             angle = 90, vjust = -0.25, size = 3) + 
    annotate("text", x = 1.15, y = 0.5, 
             label = TeX(paste0("m2B = ", params["m2B"])), 
             angle = -90, vjust = 1.5, size = 3) + 
    annotate("text", x = 0.75, y = 0.75, 
             label = TeX(paste0("m2A = ", params["m2A"])), 
             angle = 45, vjust = -0.25, size = 3) + 
    annotate("text", x = 0.25, y = 0.75, 
             label = TeX(paste0("m1B = ", params["m1B"])), 
             angle = -45, vjust = -0.25, size = 3) + 
    xlim(c(-0.4, 1.4)) + 
    coord_cartesian(ylim = c(-0.25, 1.15), clip = "off") + 
    theme_void() +
    theme(legend.position = "none",
          plot.caption = element_text(hjust = 0.5, size = 10))
  return(params_plot)
}

# Define a function for simulating the dynamics of the 
# PSF model with deSolve
psf_model <- function(time, init, params) {
  with (as.list(c(time, init, params)), {
    # description of parameters (see Bever et al. 1997)
    # N1 and N2: abundance of the plant species 1 and 2
    # p1: frequency of plant species 1; p2 = 1-pA
    # m1A, m1B, m2A, m2B: conspecific and heterospecific effects of microbial community A or B on the growth of plant 1 or 2
    # pA: frequency of the soil microbial community A
    # v: influence of plant species 2 on the microbial community relative to that of plant 1
    
    # Differential equations
    dN1 <- (m1A*pA + m1B*(1-pA))*N1
    dN2 <- (m2A*pA + m2B*(1-pA))*N2
    dp1 <- p1*(1-p1)*((m1A-m2A)*pA + (m1B-m2B)*(1-pA))
    dpA <- pA*(1-pA)*(p1-v*(1-p1))
    
    
    # Return dN1 and dN2
    return(list(c(dN1, dN2, dp1, dpA)))
  })
}

# Fig 1 --------

# Define a function for making Fig 1 so that all parameters stay cleanly in that environment
make_figure_1 <- function() {
  # Note that here, m1A technically equals the 'net' growth rate 
  # of plant 1 in microbial community A (i.e. m10+m1A)
  params <- c(m10 = 0.16, m20 = 0.16,
              m1A = 0.11, m1B = 0.26, 
              m2A = 0.27, m2B = 0.13, v = 1)
  
  time <- seq(0,50,0.1)
  
  # Run three simulations of plant growth with microbes:
  # out_pA_0 is the exponential growth of plants 1 and 2 in soil B
  init_pA_0 <- c(N1 = 3, N2 = 3, p1 = 0.5, pA = 0)
  out_pA_0 <- ode(y = init_pA_0, times = time, func = psf_model, parms = params) %>% data.frame()
  # out_pA_0 is the exponential growth of plants 1 and 2 in soil A
  init_pA_1 <- c(N1 = 3, N2 = 3, p1 = 0.5, pA = 1)
  # out_pA_05 is the growth of plants 1 and 2 in a dynamic soil
  out_pA_1 <- ode(y = init_pA_1, times = time, func = psf_model, parms = params) %>% data.frame()
  init_pA_05 <- c(N1 = 3, N2 = 7, p1 = 0.3, pA = 0.3)
  out_pA_05 <- ode(y = init_pA_05, times = time, func = psf_model, parms = params) %>% data.frame()
  # Growth of plants if there is no soil effects, calculated as
  # the exponential growth Nt = N0*e^(rt)
  no_micr <- 
    tibble(time= out_pA_05$time) %>% 
    mutate(N1 = init_pA_0["N1"]*exp(params["m10"]*time),
           N2 = init_pA_0["N2"]*exp(params["m20"]*time),
           p1 = 0,pA = 0)
  
  
  # make plot of parameters (panel A)
  # this step goes from the "composite" parameters above, to the 'actual' values
  # of miX as modifications of m10:
  params_to_plot <- c(m1A = unname(params["m1A"]-params["m10"]),
                      m1B = unname(params["m1B"]-params["m10"]),
                      m2A = unname(params["m2A"]-params["m20"]),
                      m2B = unname(params["m2B"]-params["m20"]), v=1)
  param_plot <- 
    make_params_plot(params_to_plot, scale = 6) + 
    annotate("text", x = 1.375, y = 0.5, 
             label = paste0("v = ", params["v"]), 
             angle = -90, vjust = 1.5, size = 3)  
  
  # Start making Panel B:
  # Combine simulations together for cases in which microbial
  # community is entirely A, entirely B, or asent
  for_extreme_plot <- 
    bind_rows(out_pA_1, out_pA_0, no_micr, .id = "which") %>% 
    mutate(initial_pA = case_when(which == 1 ~ "pa = 1", 
                                  which == 2 ~ "pb = 1",
                                  which == 3 ~ "none")) %>% 
    pivot_longer(N1:N2, values_to = "N", names_to = "which_sp") %>% 
    mutate(which_sp = ifelse(which_sp == "N1", "Plant 1", "Plant 2")) %>% 
    filter(time < 10.1) 
  
  # Make the panel that shows plant growth in 'extreme' soils
  panel_extreme <- 
    for_extreme_plot %>% 
    ggplot(aes(x = time, y = N, color = which_sp)) +
    geom_line(aes(linetype = initial_pA, linewidth = initial_pA)) + 
    scale_color_manual(values = c("#9970ab", "#5aae61"), labels = c("A","B"),
                       name = "species", guide = 'none') +
    scale_linewidth_manual(values = c(0.25,1.2,1.2), guide = 'none') +
    scale_linetype_manual(values = c(1,2,3),
                          labels = c("unconditioned", TeX("$S_A = 1$"),TeX("$S_B = 1$")),
                          name = "Soil community") +
    facet_wrap(.~which_sp, scales = "free") +
    guides(linetype = guide_legend(override.aes = list(linewidth = 0.5))) + 
    ylab("Abundance") + 
    scale_x_continuous(breaks = c(0,5,10)) + 
    scale_y_continuous(limits = c(0,40)) +
    theme_classic()  + 
    theme(legend.text.align = 0, axis.title = element_text(size = 10))
  
  # Panel C: plot for abundances of both plants when growing 
  # in dynamic soils
  panel_abund <- 
    out_pA_05 %>% 
    as_tibble() %>% 
    select(time, N1, N2) %>% 
    pivot_longer(N1:N2) %>%
    ggplot(aes(x = time, y = value, color = name)) +
    geom_line(linewidth = 0.9) + 
    scale_color_manual(values = c("#9970ab", "#5aae61"),
                       name = "Plant species", label = c("Plant 1", "Plant 2"), guide = "none") +
    scale_y_continuous(breaks = c(2e4, 4e4, 6e4, 8e4), labels = scales::scientific) +
    theme_classic() +
    ylab("Abundance") +
    theme(axis.title = element_text(size = 10))
  
  # Panel D: plot for frequencies of both plants when growing 
  # in dynamic soils
  panel_freq <- 
    out_pA_05 %>% 
    as_tibble() %>% 
    mutate(p2 = 1-p1) %>% 
    select(time, p1, p2) %>% 
    pivot_longer(p1:p2) %>%
    ggplot(aes(x = time, y = value, color = name)) +
    geom_line(linewidth = 1.2) + 
    scale_color_manual(values = c("#9970ab", "#5aae61"),
                       name = "Plant species", label = c("Plant 1", "Plant 2")) + 
    scale_y_continuous(limits = c(0,1), breaks = c(0, 0.5, 1)) + 
    theme_classic() +
    ylab("Frequency") +
    theme(axis.title = element_text(size = 10))
  
  # Combine panels together
  fig <- 
    {param_plot} + 
    {{panel_extreme} / 
        {{panel_abund} + 
            {panel_freq} + plot_layout(guides = "collect")}} + 
    plot_layout(widths = c(3.75,4)) + 
    plot_annotation(tag_levels = "A") & 
    theme(plot.tag = element_text(face = "bold", size = 10),
          plot.margin = unit(c(0,0,0,0), "cm"),
          legend.title = element_text(size = 10))
  
  return(fig)
}
fig1 <- make_figure_1()
ggsave("figures/fig1.pdf", plot = fig1, width = 7.25, height = 3.0, units = "in", dpi = 600)


# Fig 2 -------------
# as above, define a function for making the top panels
make_figure_2AtoD <- function() {
  params <- c(m10 = 0.16, m20 = 0.16,
              m1A = 0.02, m2A = 0.33, 
              m1B = 0.18, m2B = 0.20, v = 1)
  time <- seq(0,200,0.1)
  
  init_pA_0 <- c(N1 = 3, N2 = 3, p1 = 0.5, pA = 0)
  out_pA_0 <- ode(y = init_pA_0, times = time, func = psf_model, parms = params) %>% data.frame()
  init_pA_1 <- c(N1 = 3, N2 = 3, p1 = 0.5, pA = 1)
  out_pA_1 <- ode(y = init_pA_1, times = time, func = psf_model, parms = params) %>% data.frame()
  init_pA_05 <- c(N1 = 7, N2 = 3, p1 = 0.7, pA = 0.7)
  out_pA_05 <- ode(y = init_pA_05, times = time, func = psf_model, parms = params) %>% data.frame()
  
  
  no_micr <- 
    tibble(time= out_pA_05$time) %>% 
    mutate(N1 = init_pA_0["N1"]*exp(params["m10"]*time),
           N2 = init_pA_0["N2"]*exp(params["m20"]*time),
           p1 = 0,pA = 0)
  # Panel A 
  params_to_plot <- c(m1A = unname(params["m1A"]-params["m10"]),
                      m1B = unname(params["m1B"]-params["m10"]),
                      m2A = unname(params["m2A"]-params["m20"]),
                      m2B = unname(params["m2B"]-params["m20"]), v=1)
  
  param_plot <- make_params_plot(params_to_plot, scale = 6)
  
  # Make Panel B
  for_extreme_plot <- 
    bind_rows(out_pA_1, out_pA_0, no_micr, .id = "which") %>% 
    mutate(initial_pA = case_when(which == 1 ~ "pa = 1", 
                                  which == 2 ~ "pb = 1",
                                  which == 3 ~ "none")) %>% 
    pivot_longer(N1:N2, values_to = "N", names_to = "which_sp") %>% 
    mutate(which_sp = ifelse(which_sp == "N1", "Plant 1", "Plant 2")) %>% 
    filter(time < 10.1) 
  
  panel_extreme <- 
    for_extreme_plot %>% 
    ggplot(aes(x = time, y = N, color = which_sp)) +
    geom_line(aes(linetype = initial_pA, linewidth = initial_pA)) + 
    scale_color_manual(values = c("#9970ab", "#5aae61"), labels = c("A","B"),
                       name = "species", guide = 'none') +
    scale_linewidth_manual(values = c(0.25,1.2,1.2), guide = 'none') +
    scale_linetype_manual(values = c(1,2,3),
                          labels = c("unconditioned", TeX("$S_A = 1$"),TeX("$S_B = 1$")),
                          name = "Soil community" ) +
    facet_wrap(.~which_sp, scales = "free") +
    guides(linetype = guide_legend(override.aes = list(linewidth = 0.5))) + 
    ylab("Abundance") + 
    scale_x_continuous(breaks = c(0, 4, 8), limits = c(0,8)) +
    scale_y_continuous(limits = c(0,40)) + 
    theme_classic()  + 
    theme(legend.text.align = 0, 
          axis.title = element_text(size = 10))
  
  # Make panel C
  panel_abund <- 
    out_pA_05 %>% 
    as_tibble() %>% 
    select(time, N1, N2) %>% 
    pivot_longer(N1:N2) %>%
    ggplot(aes(x = time, y = value, color = name)) +
    geom_line(linewidth = 1.2) + 
    scale_color_manual(values = c("#9970ab", "#5aae61"),
                       name = "Plant species", label = c("Plant 1", "Plant 2")) +
    # scale_y_log10() +
    scale_y_continuous(breaks = c(0,8e17,1.6e18), labels = c("0", "7e17","1.4e18")) +
    scale_x_continuous(breaks = c(0,75,150)) +
    theme_classic() +
    ylab("Abundance") +
    theme(axis.title = element_text(size = 10))
  
  # Make panel D
  panel_freq <- 
    out_pA_05 %>% 
    as_tibble() %>% 
    mutate(p2 = 1-p1) %>% 
    select(time, p1, p2) %>% 
    pivot_longer(p1:p2) %>%
    ggplot(aes(x = time, y = value, color = name)) +
    geom_line(linewidth = 1.2) + 
    scale_color_manual(values = c("#9970ab", "#5aae61"),
                       name = "Plant species", label = c("Plant 1", "Plant 2")) + 
    scale_y_continuous(limits = c(0,1.05), breaks = c(0, 0.5, 1), expand = c(0,0)) + 
    scale_x_continuous(breaks = c(0,75,150)) +
    theme_classic() +
    ylab("Frequency") + 
    theme(axis.title = element_text(size = 10))
  
  # Combine together
  fig <- 
    {param_plot} + 
    {{panel_extreme} / 
        {{panel_abund} + 
            {panel_freq} + plot_layout(guides = "collect")}} + 
    plot_layout(widths = c(4,4)) + 
    plot_annotation(tag_levels = "A") & 
    theme(plot.tag = element_text(face = "bold", size = 10),
          plot.margin = unit(c(0,0,0,0), "cm"),
          legend.title = element_text(size = 10))
  
  return(fig)
}
fig2AtoD <- make_figure_2AtoD()

# Making Panel E is more involved, as it requires us to
# simulate the model to equilibrium and evaluate the conditions:

simulate_psf_outcomes <- function(nreps = 2000) {
  # set a seed so that we can reproduce the results exactly
  set.seed(12345)
  
  # Define a new version of the model to generate frequency-based dynamics only
  psf_model_freqOnly <- function(time, init, params) {
    with (as.list(c(time, init, params)), {
      # description of parameters (see Bever et al. 1997)
      # pA: frequency of plant species 1; pB = 1-pA
      # m1A, m1B, m2A, m2B: conspecific and heterospecific effects of microbial community A or B on the growth of plant 1 or 2
      # pAlpha: frequency of the soil microbial community A
      # v: influence of plant species 2 on the microbial community relative to that of plant 1
      
      # Differential equations
      dp1 <- p1*(1-p1)*((m1A-m2A)*pA + (m1B-m2B)*(1-pA))
      dpA <- pA*(1-pA)*(p1-v*(1-p1))
      
      
      # Return dN1 and dN2
      return(list(c(dp1, dpA)))
    })
  }
  
  
  # In this tibble, we are going to first define the parameter vectors
  # from random draws, then use each parameter vectors to run
  # the simulation. Each simulation gets run twice: once with
  # plant 1 being more frequent initially, and once with plant 2 being more frequent
  # initially. This is important for evaluating priority effects. 
  # Each substantial line of code is preceeded by a comment below:
  sim_df <-
    # Generate a vector of mIX terms, drawn from a random uniform distribution
    # from -0.5 to 0.5:
    matrix(runif(4*nreps,-0.5,0.5), nrow = nreps, ncol = 4, 
           dimnames = list(NULL, c("m1A", "m1B", "m2A", "m2B"))) %>% 
    tibble::as_tibble() %>% 
    # set v = 1 for all simulations
    mutate(v = 1) %>% 
    # Combine v and the parameters into one vector ('params_vec'),
    # and then use the resulting params_vec to run the psf model to equilibrium
    # using rootSolve::runstead(). 
    # Given the randomness of the parameters, we sometimes get combinations
    # that cannot be solved, so a call to tryCatch() ensures that the code
    # doesn't break on unstable outcomes; we just get a notice that an error occurred
    mutate(params_vec = pmap(., c),
           # root_1 runs the model with plant 2 as the more initially-frequent species (p1=0.25),
           root_1 = map(params_vec,
                        ~tryCatch(
                          {runsteady(y = c(p1=0.25, pA = 0.25), 
                                     time = c(0,Inf),
                                     func = psf_model_freqOnly,
                                     parms = .x, stol = 1e-8, rtol = 1e-6, atol = 1e-6)},
                          error = function(e) {return("error occurred")})),
           # root_2 runs the model with plant 1 as the more initially-frequent species (p1=0.75),
           root_2 = map(params_vec,
                        ~tryCatch(
                          {runsteady(y = c(p1=0.75, pA = 0.75), 
                                     time = c(0,Inf),
                                     func = psf_model_freqOnly,
                                     parms = .x, stol = 1e-8, rtol = 1e-6, atol = 1e-6)},
                          error = function(e) {return("error occurred")}
                        ))
    )
  
  # filter and manage the resulting tibble:
  sim_df_filtered <- 
    sim_df %>% 
    # throw out any parameter where either root1 or root2 failed to simulate
    filter(root_1 !="error occurred" & root_2 != "error occurred") %>% 
    mutate(
      # add some columns to calculate stabilization and fitness differneces
      # We can use these to analytically verify that our simulation works.
      stab = -0.5*(m1A-m1B-m2A+m2B),
      fdif =  0.5*(m1A+m1B-m2A-m2B),
      # For each row, identify whether analytically we exect coex/priority effects/exclusion
      coex = case_when(stab > abs(fdif) ~ "co",
                       stab < 0 & abs(stab) > abs(fdif) ~ "pe",
                       abs(stab) < abs(fdif) ~ "ex"),
      
      # Check whether the simulation 1 and 2 ran to equilibrium (steady)
      run1_steady = map_lgl(root_1, ~attr(.x, "steady")),
      run2_steady = map_lgl(root_2, ~attr(.x, "steady")),
      # Get the final frequency of plant 1 in the simulation
      run1_finalp1 = map_dbl(root_1, ~.x$y %>% pluck("p1")),
      run2_finalp1 = map_dbl(root_2, ~.x$y %>% pluck("p1"))
    ) %>% 
    # Only keep simulation runs where both runs ran to a steady state,
    # or when neither ran to a steady state (the latter corresponds to 'coexistence'),
    # since the frequencies keep oscillating
    filter(run1_steady == run2_steady) %>% 
    # Keep only important columns
    select(m1A:m2B, stab, fdif, coex, run1_steady:run2_finalp1) %>% 
    mutate(
      # make a new column of the outcome inferred from the simulation (not analytical)
      outcome = case_when(
        # if neither run ran to steady state, it is because the frequencies
        # keep oscillating; this is 'coexistence'
        !run1_steady & !run2_steady ~ "co",
        # if both ran to steady state, and the difference between the final frequencies
        # of plant 1 is minimal (less than 9e-5), this is exclusion (same plant dominates)
        run1_steady & run2_steady &
          abs(run1_finalp1 - run2_finalp1) < 9e-5 ~ "ex",
        # if both ran to steady state, but the difference between the final frequencies
        # of plant 1 is higher than threshold (9e-5), this is priority effects (outcomes
        # depend on initial frequency)
        run1_steady & run2_steady &
          abs(run1_finalp1 - run2_finalp1) > 9e-5 ~ "pe",
      ),
      # verify that our analytical solution of coexistence matches simulation
      # (confirmatory check only)
      same_ = outcome == coex
    ) 
  sim_df_filtered
}

# Running this simulation takes a long time, so I am making it simpler here:
rerun_IS_simulation = FALSE
if(rerun_IS_simulation) { 
  IS_simulation <- simulate_psf_outcomes() 
} else {
  IS_simulation <- readRDS("IS_simulation_out.Rds")
}

# Use the simulaiton output to make a plot:
fig2_panelE <-
  IS_simulation %>% 
  ggplot() + 
  geom_density(aes(x = m1A-m2A-m1B+m2B,
                   fill = outcome,
                   y = after_stat(density * n/nrow(IS_simulation))),
               bw = 0.07, alpha = 0.5, trim = TRUE, linewidth = 0.5) + 
  geom_vline(xintercept = 0, linetype = 2) +
  ylab("Density") +
  xlab(TeX("$I_S$")) + 
  annotate("text", x = -1.1, y = 0.3, label = 'coexistence',
           color = "#666633", size = 3, fontface = "bold") + 
  annotate("text", x = 0.4, y = 0.65, label = 'exclusion',
           color = "grey25", size = 3, fontface = "bold") + 
  annotate("text", x = 1.35, y = 0.25, label = 'priority\neffects',
           color = "#994455", size = 3, fontface = "bold") + 
  scale_fill_manual(values = c("#ccbb44","#bbbbbb" ,"#ee8866")) + 
  theme_classic() + 
  theme(legend.position = "none", axis.text = element_text(size = 8), 
        axis.title = element_text(size = 10))

fig2_complete <- 
  fig2AtoD/
  {plot_spacer()}/
  {plot_spacer() + fig2_panelE + plot_layout(widths = c(1,5))} + 
  plot_layout(heights = c(4,0.25,2)) + 
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 10))
ggsave("figures/fig2.pdf", plot = fig2_complete, 
       width=7.25, height = 4.5)

# Fig S2--------
# As above, this function makes Figure S2
make_figure_S2 <- function() {
  # note that the miX parameters here are identical to those in 
  # Fig 1, but now v = 5 instead of 1
  params <- c(m10 = 0.16, m20 = 0.16,
              m1A = 0.11, m1B = 0.26, m2A = 0.27, m2B = 0.13, v = 5)
  
  time <- seq(0,100,0.1)
  
  init_pA_0 <- c(N1 = 3, N2 = 3, p1 = 0.5, pA = 0)
  out_pA_0 <- ode(y = init_pA_0, times = time, func = psf_model, parms = params) %>% data.frame()
  init_pA_1 <- c(N1 = 3, N2 = 3, p1 = 0.5, pA = 1)
  out_pA_1 <- ode(y = init_pA_1, times = time, func = psf_model, parms = params) %>% data.frame()
  init_pA_05 <- c(N1 = 3, N2 = 7, p1 = 0.3, pA = 0.3)
  out_pA_05 <- ode(y = init_pA_05, times = time, func = psf_model, parms = params) %>% data.frame()
  
  no_micr <- 
    tibble(time= out_pA_05$time) %>% 
    mutate(N1 = init_pA_0["N1"]*exp(params["m10"]*time),
           N2 = init_pA_0["N2"]*exp(params["m20"]*time),
           p1 = 0,pA = 0)
  
  # Panel A
  params_to_plot <- c(m1A = unname(params["m1A"]-params["m10"]),
                      m1B = unname(params["m1B"]-params["m10"]),
                      m2A = unname(params["m2A"]-params["m20"]),
                      m2B = unname(params["m2B"]-params["m20"]), v=5)
  
  param_plot <- make_params_plot(params_to_plot, scale = 6) + 
    annotate("text", x = 1.375, y = 0.5, 
             label = paste0("v = ", params_to_plot["v"]), 
             angle = -90, vjust = 1.5, size = 3, fontface = "bold.italic")  
  
  
  # Panel B
  for_extreme_plot <- 
    bind_rows(out_pA_1, out_pA_0, no_micr, .id = "which") %>% 
    mutate(initial_pA = case_when(which == 1 ~ "pa = 1", 
                                  which == 2 ~ "pb = 1",
                                  which == 3 ~ "none")) %>% 
    pivot_longer(N1:N2, values_to = "N", names_to = "which_sp") %>% 
    mutate(which_sp = ifelse(which_sp == "N1", "Plant 1", "Plant 2")) %>% 
    filter(time < 10.1) 
  
  
  panel_extreme <- 
    for_extreme_plot %>% 
    ggplot(aes(x = time, y = N, color = which_sp)) +
    geom_line(aes(linetype = initial_pA, linewidth = initial_pA)) + 
    scale_color_manual(values = c("#9970ab", "#5aae61"), labels = c("A","B"),
                       name = "species", guide = 'none') +
    scale_linewidth_manual(values = c(0.25,1.2,1.2), guide = 'none') +
    
    scale_linetype_manual(values = c(1,2,3),
                          labels = c("unconditioned", TeX("$S_A = 1$"),TeX("$S_B = 1$")),
                          name = "Soil community") +
    facet_wrap(.~which_sp, scales = "free") +
    guides(linetype = guide_legend(override.aes = list(linewidth = 0.5))) + 
    ylab("Abundance") + 
    scale_x_continuous(breaks = c(0,5,10)) + 
    scale_y_continuous(limits = c(0,40)) + 
    theme_classic()  + 
    theme(legend.text.align = 0)
  
  
  # Panel C  
  panel_abund <- 
    out_pA_05 %>% 
    as_tibble() %>% 
    select(time, N1, N2) %>% 
    pivot_longer(N1:N2) %>%
    ggplot(aes(x = time, y = value, color = name)) +
    geom_line(linewidth = 0.9) + 
    scale_color_manual(values = c("#9970ab", "#5aae61"),
                       name = "Plant species", label = c("Plant 1", "Plant 2"), guide = "none") +
    # scale_y_log10(breaks = c(1e1, 1e3, 1e5)) +
    scale_y_continuous(breaks = c(2e8, 6e8, 1e9)) +
    theme_classic() +
    ylab("Abundance")
  
  
  # Panel D
  panel_freq <- 
    out_pA_05 %>% 
    as_tibble() %>% 
    mutate(p2 = 1-p1) %>% 
    select(time, p1, p2) %>% 
    pivot_longer(p1:p2) %>%
    ggplot(aes(x = time, y = value, color = name)) +
    geom_line(linewidth = 1.2) + 
    scale_color_manual(values = c("#9970ab", "#5aae61"),
                       name = "Plant species", label = c("Plant 1", "Plant 2")) + 
    scale_y_continuous(limits = c(0,1.05), breaks = c(0, 0.5, 1), expand = c(0,0)) + 
    theme_classic() +
    ylab("Frequency")
  
  fig <- 
    {param_plot} + 
    {{panel_extreme} / 
        {{panel_abund} + 
            {panel_freq} + plot_layout(guides = "collect")}} + 
    plot_layout(widths = c(3.5,4)) + 
    plot_annotation(tag_levels = "A") & 
    theme(plot.tag = element_text(face = "bold", size = 10),
          plot.margin = unit(c(0,0,0,0), "cm"))
  
  return(fig)
}

figS2 <- make_figure_S2()
ggsave("figures/figS2.pdf", plot = figS2, width = 8, height = 3.25)
