data <- clinical_demog_clean
install.packages("ggplot2")
library(ggplot2)

p1 <-ggplot(clinical_demog_clean, aes(x = factor(FALLER), y = BERG)) +  # factor(kaatuja) varmistaa, että käytetään kategorista akselia
  geom_boxplot(aes(color = factor(FALLER)), fill = "lightblue", alpha = 0.5) +
  labs(title = "Kaatujat vs. Ei kaatujat BERG-arvolla", 
       x = "Kaatumistulos (0 = Ei kaatuja, 1 = Kaatuja)", 
       y = "BERG-arvo") +
  theme_minimal()
print(p1)

p2 <-ggplot(clinical_demog_clean, aes(x = factor(FALLER), y = GCS_NEUROTRAX)) +  # factor(kaatuja) varmistaa, että käytetään kategorista akselia
  geom_boxplot(aes(color = factor(FALLER)), fill = "lightblue", alpha = 0.5) +
  labs(title = "FALLER VS NON-FALLER", 
       x = "FALLER (0 = NON FALLER , FALLER )", 
       y = "GCS_NEUROTRAX") +
  theme_minimal()
print(p2)

p3 <-ggplot(clinical_demog_clean, aes(x = factor(FALLER), y = EFI_EXEC_FUNC_INDEX)) +  # factor(kaatuja) varmistaa, että käytetään kategorista akselia
  geom_boxplot(aes(color = factor(FALLER)), fill = "lightblue", alpha = 0.5) +
  labs(title = "FALLER VS NON-FALLER", 
       x = "FALLER (0 = NON FALLER , FALLER )", 
       y = "EFI") +
  theme_minimal()
print(p3)

p4<- ggplot(clinical_demog_clean, aes(x = factor(FALLER), y = YEAR_FALL)) +  # factor(kaatuja) varmistaa, että käytetään kategorista akselia
  geom_boxplot(aes(color = factor(FALLER)), fill = "lightblue", alpha = 0.5) +
  labs(title = "FALLER VS NON-FALLER", 
       x = "FALLER (0 = NON FALLER , FALLER )", 
       y = "YEAR_FALL") +
  theme_minimal()
print(p4)

p5 <-ggplot(clinical_demog_clean, aes(x = factor(FALLER), y = GDS)) +  # factor(kaatuja) varmistaa, että käytetään kategorista akselia
  geom_boxplot(aes(color = factor(FALLER)), fill = "lightblue", alpha = 0.5) +
  labs(title = "FALLER VS NON-FALLER", 
       x = "FALLER (0 = NON FALLER , FALLER )", 
       y = "GDS") +
  theme_minimal()
print(p5)

p6 <- ggplot(clinical_demog_clean, aes(x = factor(FALLER), y = ABC_TOTAL_PERCENT)) +  # factor(kaatuja) varmistaa, että käytetään kategorista akselia
  geom_boxplot(aes(color = factor(FALLER)), fill = "lightblue", alpha = 0.5) +
  labs(title = "FALLER VS NON-FALLER", 
       x = "FALLER (0 = NON FALLER , FALLER )", 
       y = "ABC") +
  theme_minimal()
print(p6)

p7 <- ggplot(clinical_demog_clean, aes(x = factor(FALLER), y = SF36)) +  # factor(kaatuja) varmistaa, että käytetään kategorista akselia
  geom_boxplot(aes(color = factor(FALLER)), fill = "lightblue", alpha = 0.5) +
  labs(title = "FALLER VS NON-FALLER", 
       x = "FALLER (0 = NON FALLER , FALLER )", 
       y = "SF36") +
  theme_minimal()
print(p7)

p8 <- ggplot(clinical_demog_clean, aes(x = factor(FALLER), y = PASE)) +  # factor(kaatuja) varmistaa, että käytetään kategorista akselia
  geom_boxplot(aes(color = factor(FALLER)), fill = "lightblue", alpha = 0.5) +
  labs(title = "FALLER VS NON-FALLER", 
       x = "FALLER (0 = NON FALLER , FALLER )", 
       y = "PASE") +
  theme_minimal()
print(p8)

p9 <- ggplot(clinical_demog_clean, aes(x = factor(FALLER), y = MMSE)) +  # factor(kaatuja) varmistaa, että käytetään kategorista akselia
  geom_boxplot(aes(color = factor(FALLER)), fill = "lightblue", alpha = 0.5) +
  labs(title = "FALLER VS NON-FALLER", 
       x = "FALLER (0 = NON FALLER , FALLER )", 
       y = "MMSE") +
  theme_minimal()
print(p9)

p10 <- ggplot(clinical_demog_clean, aes(x = factor(FALLER), y = MOCA)) +  # factor(kaatuja) varmistaa, että käytetään kategorista akselia
  geom_boxplot(aes(color = factor(FALLER)), fill = "lightblue", alpha = 0.5) +
  labs(title = "FALLER VS NON-FALLER", 
       x = "FALLER (0 = NON FALLER , FALLER )", 
       y = "MOCA") +
  theme_minimal()
print(p10)

p11 <- ggplot(clinical_demog_clean, aes(x = factor(FALLER), y = FAB)) +  # factor(kaatuja) varmistaa, että käytetään kategorista akselia
  geom_boxplot(aes(color = factor(FALLER)), fill = "lightblue", alpha = 0.5) +
  labs(title = "FALLER VS NON-FALLER", 
       x = "FALLER (0 = NON FALLER , FALLER )", 
       y = "FAB") +
  theme_minimal()
print(p11)

p12<- ggplot(clinical_demog_clean, aes(x = factor(FALLER), y = TMT_A)) +  # factor(kaatuja) varmistaa, että käytetään kategorista akselia
  geom_boxplot(aes(color = factor(FALLER)), fill = "lightblue", alpha = 0.5) +
  labs(title = "FALLER VS NON-FALLER", 
       x = "FALLER (0 = NON FALLER , FALLER )", 
       y = "TMT_A") +
  theme_minimal()
print(p12)

p13 <- ggplot(clinical_demog_clean, aes(x = factor(FALLER), y = TMT_B)) +  # factor(kaatuja) varmistaa, että käytetään kategorista akselia
  geom_boxplot(aes(color = factor(FALLER)), fill = "lightblue", alpha = 0.5) +
  labs(title = "FALLER VS NON-FALLER", 
       x = "FALLER (0 = NON FALLER , FALLER )", 
       y = "TMT_B") +
  theme_minimal()
print(p13)

p14 <- ggplot(clinical_demog_clean, aes(x = factor(FALLER), y = TUG)) +  # factor(kaatuja) varmistaa, että käytetään kategorista akselia
  geom_boxplot(aes(color = factor(FALLER)), fill = "lightblue", alpha = 0.5) +
  labs(title = "FALLER VS NON-FALLER", 
       x = "FALLER (0 = NON FALLER , FALLER )", 
       y = "TUG") +
  theme_minimal()
print(p14)

p15 <- ggplot(clinical_demog_clean, aes(x = factor(FALLER), y = FSST)) +  # factor(kaatuja) varmistaa, että käytetään kategorista akselia
  geom_boxplot(aes(color = factor(FALLER)), fill = "lightblue", alpha = 0.5) +
  labs(title = "FALLER VS NON-FALLER", 
       x = "FALLER (0 = NON FALLER ,1= FALLER )", 
       y = "FSST") +
  theme_minimal()
print(p15)

p16 <-ggplot(clinical_demog_clean, aes(x = factor(FALLER), y = BERG)) +  # factor(kaatuja) varmistaa, että käytetään kategorista akselia
  geom_boxplot(aes(color = factor(FALLER)), fill = "lightblue", alpha = 0.5) +
  labs(title = "FALLER VS NON-FALLER", 
       x = "FALLER (0 = NON FALLER ,1= FALLER )", 
       y = "BERG") +
  theme_minimal()
print(p16)

p17<- ggplot(clinical_demog_clean, aes(x = factor(FALLER), y = DGI)) +  # factor(kaatuja) varmistaa, että käytetään kategorista akselia
  geom_boxplot(aes(color = factor(FALLER)), fill = "lightblue", alpha = 0.5) +
  labs(title = "FALLER VS NON-FALLER", 
       x = "FALLER (0 = NON FALLER ,1= FALLER )", 
       y = "DGI") +
  theme_minimal()
print(p17)

p18 <- ggplot(clinical_demog_clean, aes(x = factor(FALLER), y = DGI_STAIRS)) +  # factor(kaatuja) varmistaa, että käytetään kategorista akselia
  geom_boxplot(aes(color = factor(FALLER)), fill = "lightblue", alpha = 0.5) +
  labs(title = "FALLER VS NON-FALLER", 
       x = "FALLER (0 = NON FALLER ,1= FALLER )", 
       y = "DGI_STAIRS") +
  theme_minimal()
print(p18)

p19 <- ggplot(clinical_demog_clean, aes(x = factor(FALLER), y = BASE_VELOCITY)) +  # factor(kaatuja) varmistaa, että käytetään kategorista akselia
  geom_boxplot(aes(color = factor(FALLER)), fill = "lightblue", alpha = 0.5) +
  labs(title = "FALLER VS NON-FALLER", 
       x = "FALLER (0 = NON FALLER ,1= FALLER )", 
       y = "BASE_VELOCITY") +
  theme_minimal()
print(p19)

p20 <- ggplot(clinical_demog_clean, aes(x = factor(FALLER), y = S3_VELOCITY)) +  # factor(kaatuja) varmistaa, että käytetään kategorista akselia
  geom_boxplot(aes(color = factor(FALLER)), fill = "lightblue", alpha = 0.5) +
  labs(title = "FALLER VS NON-FALLER", 
       x = "FALLER (0 = NON FALLER ,1= FALLER )", 
       y = "S3_VELOCITY") +
  theme_minimal()
print(p20)

p21 <- ggplot(clinical_demog_clean, aes(x = factor(FALLER), y = AGE)) +  # factor(kaatuja) varmistaa, että käytetään kategorista akselia
  geom_boxplot(aes(color = factor(FALLER)), fill = "lightblue", alpha = 0.5) +
  labs(title = "FALLER VS NON-FALLER", 
       x = "FALLER (0 = NON FALLER ,1= FALLER )", 
       y = "AGE") +
  theme_minimal()
print(p21)


## Keskiarvo isompi kaatujilla: YEAR_FALL, GDS

## Keskiarvo isompi ei- kaatujilla: GCS_Neurotrax, SF36, MMSE(vinojakauma), MOCA, FAB
## FFST (vinojakauma?), BERG, DGI(vinojakauma?), BASE_VELOCITY, S3_VELOCITY

## Pieni ero keskiarvoilla: EFI, ABC (vaikea verrata), PASE, 
##TMT_A ja TMT_B (kaatujilla hieman isompi ka), 
##TUG(VAIKEA VERRATA), DGI_Stairs(vaikea verrata), AGE