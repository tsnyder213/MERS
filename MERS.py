import matplotlib.pyplot as plt
import numpy as np

P = 170  # P = recruitment rate of susceptible individuals
beta = 0.122  # beta = effective contact rate
subE = 0.2996  # subE = reduction factor in transmission rate by exposed per day
subJ = 0.0899  # subJ = reduction factor in transmission rate by isolated per day
k = 0.1026  # k = rate of development of clinical symptoms in exposed population
mu = 0.000034  # mu = natural death rate
d1 = 0.0294  # d1 = disease-related death rate of I (infected)
d2 = 0.0227  # d2 = disease-related death rate of J (isolated)
sig1 = 0.0433  # sig1 = recovery rate of I (infected)
sig2 = 0.0475  # sig2 = recovery rate of J (isolated)
# gamma = isolation rate (control variable)
t = 1500


def interaction(S, E, I, J, R):
    return S*(((beta*I) + (subE*beta*E) + (subJ*beta*J))/(S + E + I + J + R))


def Iplus(S, E, I, J, R, gamma):
    return I + ((k*E) - (gamma*I + sig1*I + d1*I + mu*I))


def Splus(S, E, I, J, R):
    return S + (P - (interaction(S, E, I, J, R)) - (mu*S))


def Eplus(S, E, I, J, R):
    return E + (interaction(S, E, I, J, R) - (k+mu)*E)


def Jplus(S, E, I, J, R, gamma):
    return J + (gamma*I - (sig2*J + d2*J + mu*J))


def Rplus(S, E, I, J, R):
    return R + ((sig1*I) + (sig2*J) - (mu*R))


def dp(gamma, S0=4999990, E0=0, I0=10, J0=0, R0=0):
    I = [I0]
    S = [S0]
    E = [E0]
    R = [R0]
    J = [J0]
    deaths = [0]
    for i in range(1, t+1):
        deaths.append(d1*I[i-1] + d2*J[i-1])
        R.append(Rplus(S[i - 1], E[i - 1], I[i - 1], J[i - 1], R[i - 1]))
        I.append(Iplus(S[i - 1], E[i - 1], I[i - 1], J[i - 1], R[i - 1], gamma))
        S.append(Splus(S[i - 1], E[i - 1], I[i - 1], J[i - 1], R[i - 1]))
        E.append(Eplus(S[i - 1], E[i - 1], I[i - 1], J[i - 1], R[i - 1]))
        J.append(Jplus(S[i - 1], E[i - 1], I[i - 1], J[i - 1], R[i - 1], gamma))

    return [S,E,I,J,R,deaths]


# Run the simulation for the 3 gamma values
[S, E, I, J, R, deaths] = dp(gamma=0.1501)
[S2, E2, I2, J2, R2, deaths2] = dp(gamma=0.16)
[S3, E3, I3, J3, R3, deaths3] = dp(gamma=0.17)

# Figure 1 = Infections
plt.figure()

plt.title("MERS Infections")
plt.ylabel('Number of Infected Individuals')
plt.xlabel('Time Steps (days)')

xpoints = np.array(range(0,1501))
ypoints = np.array(np.cumsum(I))
y2points = np.array(np.cumsum(I2))
y3points = np.array(np.cumsum(I3))

plt.plot(xpoints, ypoints, label="gamma = 0.1501")
plt.plot(xpoints, y2points, label="gamma = 0.16")
plt.plot(xpoints, y3points, label="gamma = 0.17")
plt.legend()

# Figure 2 = Recoveries
plt.figure()

plt.title("MERS Recoveries")
plt.ylabel('Number of Recovered Individuals')
plt.xlabel('Time Steps (days)')

xpoints = np.array(range(0,1501))
ypoints = np.array(R)
y2points = np.array(R2)
y3points = np.array(R3)

plt.plot(xpoints, ypoints, label="gamma = 0.1501")
plt.plot(xpoints, y2points, label="gamma = 0.16")
plt.plot(xpoints, y3points, label="gamma = 0.17")
plt.legend()

# Figure 3 = Deaths
plt.figure()

plt.title("MERS Deaths")
plt.ylabel('Number of Disease-Related Deaths')
plt.xlabel('Time Steps (days)')

xpoints = np.array(range(0,1501))
ypoints = np.array(np.cumsum(deaths))
y2points = np.array(np.cumsum(deaths2))
y3points = np.array(np.cumsum(deaths3))

plt.plot(xpoints, ypoints, label="gamma = 0.1501")
plt.plot(xpoints, y2points, label="gamma = 0.16")
plt.plot(xpoints, y3points, label="gamma = 0.17")
plt.legend()

plt.show()




