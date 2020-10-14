```@meta
CurrentModule = Uchiyama
```

# Uchiyama

Uchiyama’s particle motion is computed using an Event-Driven algorithm. The Event-Driven method is concerned with the times at which events, in this case collisions, take place. The algorithm is the following: you establish a list of the collisions to come if the particles moved only in straight lines. This list is ordered according to the time at which the collisions will take place, the first element being the closest collision in time. The particles are then displaced until this time and the collision is carried out. The list of future collisions is then updated. Indeed, some of them may be invalidated by the collision that has just taken place since the velocities and therefore the direction of the two particles involved have changed, new ones may also become possible for the same reasons. Then, again, we move to the nearest collision in time and so on ...

```@index
```

```@autodocs
Modules = [Uchiyama]
```

